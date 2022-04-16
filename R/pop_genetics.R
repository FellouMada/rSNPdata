#' Calculate Fst
#' @param snpdata SNPdata object
#' @param groups a vector of string. the differentiation will be estimated between these groups of samples
#' @param from the metadata column from which these groups belong to
#' @return SNPdata object with an extra field: Fst
#' @usage  calculate_wcFst(snpdata, groups=c("Senegal","Gambia"), from="Country")
#' @export
calculate_wcFst = function(snpdata, groups=NULL, from=NULL){
    if(is.null(groups) & is.null(from)){
        stop("Please provide a vector of groups to be compared and the metadata column of interest")
    }else if(is.null(groups) & !is.null(from)){
        groups = as.character(unique(snpdata$meta[[from]]))
    }else if(!is.null(groups) & !is.null(from) & (any(!(groups %in% unique(snpdata$meta[[from]]))))){
        stop("not all specified groups belong to the ", from, "column of the metadata table")
    }
    system(sprintf("tabix %s", snpdata$vcf))
    snpdata$Fst=list()
    for(i in 1:(length(groups)-1)){
        idx1 = which(snpdata$meta[[from]]==groups[i])
        idx1 = paste(idx1, collapse = ",")
        for(j in (i+1):length(groups)){
            idx2 = which(snpdata$meta[[from]]==groups[j])
            idx2 = paste(idx2, collapse = ",")
            out = paste0(dirname(snpdata$vcf),"/out.wc.fst")
            pout = paste0(dirname(snpdata$vcf),"/out.wc.fst.pvalues")
            system(sprintf("wcFst --target %s --background %s --file %s --deltaaf 0 --type GT > %s", idx1, idx2, snpdata$vcf, out))
            system(sprintf("pFst --target %s --background %s --file %s --deltaaf 0 --type GT > %s", idx1, idx2, snpdata$vcf, pout))
            out = fread(out)
            pout = fread(pout)
            setkeyv(out,c("V1","V2")); setkeyv(pout,c("V1","V2"))
            tmp = out[pout, nomatch=0]
            names(tmp) = c("Chrom","Pos","AF_in_target","AF_in_background","wcFst","wcFst_pvalue")
            idx = which(tmp$wcFst<0)
            if(length(idx)>0){
                tmp$wcFst[idx]=0
            }
            tmp$wcFst_Adj_pvalue_BH = p.adjust(tmp$wcFst_pvalue, method="BH")
            snpdata$Fst[[paste0(groups[i],"_vs_",groups[j])]] = tmp
        }
    }
    snpdata
}

#' Calculate LD R^2 between pairs of loci
#' @param snpdata SNPdata object
#' @param min.r2 the minimum r2 value below which the LD value is not reported
#' @param inter.chrom whether to calculate inter-chromosomal LD. FALSE by default
#' @param chroms a vector of chromosome names for which LD should be calculated
#' @return SNPdata object with an extra field: LD
#' @usage  calculate_LD(snpdata, min.r2=0.2, inter.chrom=FALSE, chroms=c("Pf3D7_04_v3","Pf3D7_05_v3"))
#' @details the output file from LD calculation could be large. In order to reduce the size of that file, it's recommended to specify the list of chromosomes for which LD should be calculated using the `chroms` option
#' @export
calculate_LD = function(snpdata, min.r2=0.2, inter.chrom=FALSE, chroms=NULL){
    out = paste0(dirname(snpdata$vcf),"/tmp_ld")
    if(inter.chrom){
        cat("inter-chromosomal LD will be calculated between sites on ",paste(chroms, collapse = ","))
        system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --interchrom-geno-r2", snpdata$vcf, out, min.r2))
    }else{
        system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --geno-r2", snpdata$vcf, out, min.r2))
    }
    out = paste0(out,".geno.ld")
    system(sprintf("bgzip %s", out))
    ld = fread(out, nThread = 4)
    if(!is.null(chroms)){
        idx1 = which(ld$CHR1 %in% chroms)
        idx2 = which(ld$CHR2 %in% chroms)
        idx = unique(c(idx1, idx2))
        ld = l[idx,]
    }

    snpdata$LD=ld
}

#' Calculate IBS matrix between all pairs of isolates
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="GT"
#' @return SNPdata object with an extra field: IBS
#' @usage  get_ibs(snpdata, mat.name="GT")
#' @export
get_ibs = function(snpdata, mat.name="GT"){
    if(missing(snpdata[[mat.name]])){
        stop("specified genotype matrix does not exist")
    }
    X = snpdata[[mat.name]]
    y = matrix(NA, nrow = nrow(X), ncol = nrow(X))
    colnames(y) = colnames(X)
    for(i in 1:nrow(X)){
        for(j in 1:nrow(X)){
            m = X[i,]-X[j,]
            y[i,j] = 1-(length(which(m==0))/ncol(X))
        }
    }
    snpdata$IBS=y[upper.tri(y)]=NA
    snpdata
}

get_ibd = function(snpdata){

}

#' Calculate iR index to detect loci with excess of IBD sharing
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="Phased"
#' @param family the name of the column, in the metadata table, to be used to represent the sample's population
#' @return SNPdata object with an extra field: iR
#' @usage  calculate_iR(snpdata, mat.name="Phased", family="Location")
#' @export
calculate_iR = function(snpdata, mat.name="Phased", family="Location"){
    if(!(family %in% names(snpdata$meta))){
        stop("No column name ",family," in the metadata table")
    }
    if(mat.name=="GT"){
        cat("phasing the mixed genotypes\n")
        snpdata = phase_mixed_genotypes(snpdata, nsim=100)
    }

}

make_ped = function(mat, metadata, family){
    mat[mat==1]=2
    mat[mat==0]=1
    mat[is.na(mat)]=0
    mat=t(mat)

    new.genotype = matrix(-9 , nrow = dim(mat)[1], ncol = (dim(mat)[2]*2) )
    j=1
    k=1
    while(j<=ncol(mat)){
        new.genotype[,k] = mat[,j]
        new.genotype[,(k+1)] = mat[,j]
        j=j+1
        k=k+2
    }

    PED6 = data.frame(cbind(metadata[[family]],phenotype[[sample]]), FatherID = 0, MotherID = 0, Sex = 0, Phenotype = -9, stringsAsFactors = FALSE)
    ped = data.frame(cbind(PED6, data.frame(new.genotype)))
    ped
}

make_map = function(details){
    c1 = details$Chrom
    c4 = details$Pos
    get_xme = function(x){as.character(unlist(strsplit(x,"_"))[2])}
    xme = as.numeric(as.character(lapply(c1, get_xme)))
    c2 = paste0(LETTERS[xme], c4)
    c3 = c4/12000
    map = data.frame(chrom = c1, post = c2, gd = c3, pos = c4, stringsAsFactors = FALSE)
    map
}


