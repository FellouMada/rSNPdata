#' Calculate Weir & Cockerham's Fst
#' @param snpdata SNPdata object
#' @param groups a vector of string. the differentiation will be estimated between these groups of samples
#' @param from the metadata column from which these groups belong to
#' @return SNPdata object with an extra field: Fst. This is a list of data frames that contain the result from a specific comparison. Every data frame contains N rows (where N is the number of loci) and the following 7 columns:
#' \enumerate{
#' \item the chromosome ID (type \code{"character"})
#' \item the SNPs positions (type \code{"numeric"})
#' \item the allele frequency in the first group (type \code{"numeric"})
#' \item the allele frequency in the second group (type \code{"numeric"})
#' \item the resulting Fst values (type \code{"numeric"})
#' \item the p-values associated with the Fst results (type \code{"numeric"})
#' \item the p-values corrected for multiple testing using the Benjamini-Hochberg method (type \code{"numeric"})
#' }
#' @usage  calculate_wcFst(snpdata, groups=c("Senegal","Gambia"), from="Country")
#' @export
calculate_wcFst = function(snpdata, from=NULL, groups=NULL){
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
calculate_IBS = function(snpdata, mat.name="GT"){
    if(!(mat.name %in% names(snpdata))){
        stop("specified genotype matrix does not exist")
    }
    X = t(snpdata[[mat.name]])
    y = matrix(NA, nrow = nrow(X), ncol = nrow(X))
    # colnames(y) = rownames(y) = rownames(X)
    pb = txtProgressBar(min = 0, max = nrow(X), initial = 0,style = 3)
    for(i in 1:nrow(X)){
        for(j in 1:nrow(X)){
            m = X[i,]-X[j,]
            y[i,j] = 1-(length(which(m==0))/ncol(X))
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    y[upper.tri(y)]=NA
    colnames(y) = rownames(X); rownames(y) = rownames(X)
    snpdata$IBS=y
    snpdata
}

calculate_IBD = function(snpdata){

}

#' Calculate iR index to detect loci with excess of IBD sharing
#' @param snpdata SNPdata object
#' @param mat.name the name of the genotype table to be used. default="Phased"
#' @param family the name of the column, in the metadata table, to be used to represent the sample's population
#' @param number.cores the number of cores to be used. default=4
#' @return SNPdata object with an extra field: iR
#' @usage  calculate_iR(snpdata, mat.name="Phased", family="Location", number.cores=4)
#' @export
calculate_iR = function(snpdata, mat.name="Phased", family="Location", number.cores=4){
    if(!(family %in% names(snpdata$meta))){
        stop("No column name ",family," in the metadata table")
    }
    if(mat.name=="GT"){
        cat("phasing the mixed genotypes\n")
        snpdata = phase_mixed_genotypes(snpdata, nsim=100)
    }
    ped = make_ped(snpdata[[mat.name]], snpdata$meta, family)
    ped$sex=as.numeric(ped$sex)
    map = make_map(snpdata$details)
    ped.map = list(ped,map)
    my.geno = getGenotypes(ped.map, reference.ped.map=NULL, maf=0.01, isolate.max.missing=0.2, snp.max.missing=0.2, chromosomes=NULL, input.map.distance="cM", reference.map.distance="cM")
    my.param = getIBDparameters(ped.genotypes = my.geno, number.cores = number.cores)
    my.ibd = getIBDsegments(ped.genotypes = my.geno,parameters = my.param, number.cores = number.cores, minimum.snps = 20, minimum.length.bp = 50000,error = 0.001)
    my.matrix = getIBDmatrix(ped.genotypes = my.geno, ibd.segments = my.ibd)
    my.iR = getIBDiR(ped.genotypes = my.geno,ibd.matrix = my.matrix,groups = NULL)
    pvalues = as.numeric(as.character(lapply(my.iR$log10_pvalue, get.pvalue)))
    my.iR$log10_pvalue = pvalues
    names(my.iR)[8] = "pvalue"
    my.iR$adj_pvalue_BH = p.adjust(my.iR$pvalues, method = "BH")
    if(!("iR" %in% names(snpdata))){
        snpdata$iR=list()
    }
    groups = unique(snpdata$meta[[family]])
    snpdata$iR[[paste0(groups[1],"_vs_",groups[2])]] = my.iR
    snpdata
}

get.pvalue = function(x){10^-x}

make_ped = function(mat, metadata, family){
    mat[mat==1]=2
    mat[mat==0]=1
    mat[is.na(mat)]=0
    mat=t(mat)

    new.genotype = matrix(-9 , nrow = dim(mat)[1], ncol = ((dim(mat)[2])*2) )
    j=1
    k=1
    while(j<=ncol(mat)){
        # print(j)
        new.genotype[,k] = mat[,j]
        new.genotype[,(k+1)] = mat[,j]
        j=j+1
        k=k+2
    }

    PED6 = data.frame(cbind(metadata[[family]],metadata$sample, fatherID = 0, motherID = 0, sex = metadata$COI, phenotype = -9), stringsAsFactors = FALSE)
    names(PED6)[1:2] = c("pop","sample")
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

haplotypeToGenotype = function(haplotypes, moi){
    genotypes = matrix(-9, dim(haplotypes)[1], dim(haplotypes)[2]/2)
    for(i in 1:dim(haplotypes)[1]){
        #print(i)
        j=1
        k=1
        while(j <= dim(haplotypes)[2]){
            if(haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 1) genotypes[i,k] = 0
            if(haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 2) genotypes[i,k] = 2
            if(haplotypes[i,j] == 0 && haplotypes[i,(j+1)] == 0) genotypes[i,k] = -1
            if((moi[i] == 1) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = -1
            if((moi[i] == 2) && ((haplotypes[i,j] == 1 && haplotypes[i,(j+1)] == 2) || (haplotypes[i,j] == 2 && haplotypes[i,(j+1)] == 1))) genotypes[i,k] = 1
            if(haplotypes[i,j] != 0 && haplotypes[i,j] != 1 && haplotypes[i,j] != 2) genotypes[i,k] = -1
            if(haplotypes[i,(j+1)] != 0 && haplotypes[i,(j+1)] != 1 && haplotypes[i,(j+1)] != 2) genotypes[i,k] = -1
            j=j+2
            k=k+1
        }
    }
    return(t(genotypes))
}

calculatePopAlleleFreq = function(genotypes, moi){
    pop_allele_freqs = vector('numeric',dim(genotypes)[1])
    number_isolates = dim(genotypes)[2]
    number_snps = dim(genotypes)[1]

    for (t in 1:number_snps){
        #print(t)
        A = 0
        B = 0
        for (i in 1:number_isolates)
        {
            if (genotypes[t,i] == 0 && moi[i] == 1)  A=A+1
            if (genotypes[t,i] == 0 && moi[i] == 2)  A=A+2
            if (genotypes[t,i] == 1 && moi[i] == 2){ A=A+1; B=B+1 }
            if (genotypes[t,i] == 2 && moi[i] == 1)  B=B+1
            if (genotypes[t,i] == 2 && moi[i] == 2)  B=B+1
        }
        if (A + B == 0) pop_allele_freqs[t] = -1
        else pop_allele_freqs[t] = A/(A+B)
    }
    return(pop_allele_freqs)
}

calculateMissingness = function(genotypes) {
    proportion_missing=vector('numeric',dim(genotypes)[2])
    number_snps = dim(genotypes)[2]
    number_isolates = dim(genotypes)[1]
    number_snps_1 = dim(genotypes)[1]

    for (i in 1:number_snps){
        #print(i)
        number_missing = 0.0
        for (j in 1:number_isolates) {
            if(genotypes[j,i] == -1 ) number_missing = number_missing+1;
        }
        proportion_missing[i] = number_missing/number_snps_1;
    }
    return(proportion_missing)
}

getGenotypes = function(ped.map, reference.ped.map=NULL, maf=0.01, isolate.max.missing=0.1, snp.max.missing=0.1, chromosomes=NULL, input.map.distance="cM", reference.map.distance="cM"){
    # check input PED and MAP files
    if (!is.list(ped.map) | length(ped.map) != 2) stop ("'ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
    input.ped <- as.data.frame(ped.map[[1]])
    input.map <- as.data.frame(ped.map[[2]])
    if (!is.data.frame(input.ped)) stop ("'ped.map' has incorrect format - PED is not a data.frame")
    if (!is.data.frame(input.map)) stop ("'ped.map' has incorrect format - MAP is not a data.frame")

    # check the PED and MAP files have the same number of SNPs
    if (ncol(input.ped) != (2*nrow(input.map)+6)) stop ("'ped.map' has incorrect format - PED and MAP must have the same number of SNPs")

    # check the MAP file has 4 coloumns
    if (ncol(input.map) != 4) stop ("'ped.map' has incorrect format - MAP must have four columns")
    colnames(input.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
    if (!is.numeric(input.map[,"pos_M"])) stop ("'ped.map' has incorrect format - genetic-map positions in MAP file are non-numeric")
    if (!is.numeric(input.map[,"pos_bp"])) stop ("'ped.map' has incorrect format - base-pair positions in MAP file are non-numeric")
    if (is.factor(input.map[,"chr"]))
        input.map[,"chr"] <- as.character(input.map[,"chr"])
    if (is.factor(input.map[,"snp_id"]))
        input.map[,"snp_id"] <- as.character(input.map[,"snp_id"])


    # check PED for factors
    # input.ped = as.matrix(input.ped)
    # input.map = as.matrix(input.map)
    if (is.factor(input.ped[,1]))
        input.ped[,1] <- as.character(input.ped[,1])
    if (is.factor(input.ped[,2]))
        input.ped[,2] <- as.character(input.ped[,2])


    # check reference data
    if (!is.null(reference.ped.map)) {
        if (!is.list(reference.ped.map) | length(reference.ped.map) != 2) stop ("'reference.ped.map' must be a list containing 2 objects: 'PED' and 'MAP'")
        reference.ped <- reference.ped.map[[1]]
        reference.map <- reference.ped.map[[2]]
        if (!is.data.frame(reference.ped)) stop ("'reference.ped.map' has incorrect format - PED is not a data.frame")
        if (!is.data.frame(reference.map)) stop ("'reference.ped.map' has incorrect format - MAP is not a data.frame")

        # check the PED and MAP files have the same number of SNPs
        if (ncol(reference.ped) != (2*nrow(reference.map)+6)) stop ("'reference.ped.map' has incorrect format - PED and MAP must have the same number of SNPs")

        # check the MAP file has 4 coloumns
        if (ncol(reference.map) != 4) stop ("'reference.ped.map' has incorrect format - MAP must have four columns")
        colnames(reference.map) <- c("chr", "snp_id", "pos_M", "pos_bp")
        if (!is.numeric(reference.map[,"pos_M"])) stop ("'reference.ped.map' has incorrect format - genetic-map positions in reference MAP file are non-numeric")
        if (!is.numeric(reference.map[,"pos_bp"])) stop ("'reference.ped.map' has incorrect format - base-pair positions in reference MAP file are non-numeric")
        if (is.factor(reference.map[,"chr"]))
            reference.map[,"chr"] <- as.character(reference.map[,"chr"])
        if (is.factor(reference.map[,"snp_id"]))
            reference.map[,"snp_id"] <- as.character(reference.map[,"snp_id"])

        # check PED for factors
        if (is.factor(reference.ped[,1]))
            reference.ped[,1] <- as.character(reference.ped[,1])
        if (is.factor(reference.ped[,2]))
            reference.ped[,2] <- as.character(reference.ped[,2])
    }

    # check maf
    if (!is.vector(maf)) stop ("'maf' has incorrect format - must be a vector")
    if (!is.numeric(maf)) stop ("'maf' has incorrect format - must be numeric")
    if (length(maf) != 1) stop ("'maf' has incorrect format - must be a single numeric value")
    if (maf > 1 | maf < 0) stop ("'maf' has incorrect format - must be between 0 and 1")

    # check isolate.max.missing
    if (!is.vector(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(isolate.max.missing)) stop ("'isolate.max.missing' has incorrect format - must be numeric")
    if (length(isolate.max.missing) != 1) stop ("'isolate.max.missing' has incorrect format - must be a single numeric value")
    if (isolate.max.missing > 1 | isolate.max.missing < 0) stop ("'isolate.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")

    # check snp.max.missing
    if (!is.vector(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be a vector")
    if (!is.numeric(snp.max.missing)) stop ("'snp.max.missing' has incorrect format - must be numeric")
    if (length(snp.max.missing) != 1) stop ("'snp.max.missing' has incorrect format - must be a single numeric value")
    if (snp.max.missing > 1 | snp.max.missing < 0) stop ("'snp.max.missing' has incorrect format - must be between 0 and 1 (inclusive)")

    # check chromosomes
    if (!is.null(chromosomes)) {
        if (!is.vector(chromosomes)) stop ("'chromosomes' has incorrect format - must be a vector")
        if(!all(chromosomes %in% input.map[,"chr"]))
            stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% input.map[,"chr"])]," not in 'ped.map'\n")))
        if (!is.null(reference.ped.map)) {
            if(!all(chromosomes %in% reference.map[,"chr"]))
                stop(paste0("chromosome ",paste0(chromosomes[!(chromosomes %in% reference.map[,"chr"])]," not in 'reference.ped.map'\n")))
        }
    } else
        chromosomes <- unique(as.character(input.map[,"chr"]))

    # check input map distance
    if (!is.vector(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be a vector")
    if (!is.character(input.map.distance)) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (length(input.map.distance) != 1) stop ("'input.map.distance' has incorrect format - must be either character M or cM")
    if (input.map.distance != "M" & input.map.distance != "cM")
        stop ("'input.map.distance' has incorrect format - must be either M or cM")
    if (input.map.distance == "cM") {
        input.map[,"pos_M"] <- input.map[,"pos_M"]/100
    }

    # check reference map distance
    if (!is.vector(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be a vector")
    if (!is.character(reference.map.distance)) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (length(reference.map.distance) != 1) stop ("'reference.map.distance' has incorrect format - must be either character M or cM")
    if (reference.map.distance != "M" & reference.map.distance != "cM")
        stop ("'reference.map.distance' has incorrect format - must be either M or cM")
    if (reference.map.distance == "cM" & !is.null(reference.ped.map)) {
        reference.map[,"pos_M"] <- reference.map[,"pos_M"]/100
    }

    # begin data filtering

    # create new isolate IDs from PED FIDs and IIDs
    isolate.names <- paste(input.ped[,1], input.ped[,2], sep="/")
    if (any(duplicated(isolate.names)))
        stop ("duplicate sample IDs found")

    # merge input data with reference data
    if (!is.null(reference.ped.map)) {
        input.map.v1      <- cbind(1:nrow(input.map), input.map)
        reference.map.v1  <- cbind(1:nrow(reference.map), reference.map)
        input.map.v1      <- merge(input.map.v1, reference.map.v1, by.x="snp_id", by.y="snp_id")
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map'")
        input.map.v1  <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
        if (!is.null(chromosomes))
            input.map.v1 <- input.map.v1[input.map.v1[,"chr.x"] %in% chromosomes,]
        if (nrow(input.map.v1) == 0)
            stop ("no SNPs remaining after merging 'ped.map' and 'reference.ped.map' for selected chromosomes")
        input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
        input.ped.columns <- input.ped.columns[order(input.ped.columns)]
        input.ped.v1      <- input.ped[,input.ped.columns]
        reference.ped.columns  <- c(1:6, 2*input.map.v1[,"1:nrow(reference.map)"] + 5, 2*input.map.v1[,"1:nrow(reference.map)"] + 6)
        reference.ped.columns  <- reference.ped.columns[order(reference.ped.columns)]
        reference.ped.v1       <- reference.ped[,reference.ped.columns]
        input.map.v2           <- input.map.v1[,c("chr.x", "snp_id", "pos_M.x", "pos_bp.x")]
        colnames(input.map.v2) <- c("chr", "snp_id", "pos_M","pos_bp")
    } else {
        if (!is.null(chromosomes)) {
            input.map.v1 <- cbind(1:nrow(input.map), input.map)
            input.map.v1 <- input.map.v1[input.map.v1[,"chr"] %in% chromosomes,]
            if (nrow(input.map.v1) == 0)
                stop ("no SNPs remaining after subsetting 'ped.map'by selected chromosomes")
            input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
            input.ped.columns <- input.ped.columns[order(input.ped.columns)]
            input.ped.v1      <- input.ped[,input.ped.columns]
            input.map.v2      <- input.map.v1[,c("chr", "snp_id", "pos_M", "pos_bp")]
        } else {
            input.map.v2 <- input.map
            input.ped.v1 <- input.ped
        }
    }

    # call genotypes
    input.matrix        <- as.matrix(input.ped.v1[,7:ncol(input.ped.v1)])
    input.genders       <- input.ped.v1[,5]
    input.genotypes.v0  <- cbind(input.map.v2, haplotypeToGenotype(input.matrix, input.genders))
    if (!is.null(reference.ped.map)) {
        reference.matrix       <- as.matrix(reference.ped.v1[,7:ncol(reference.ped.v1)])
        reference.genders      <- reference.ped.v1[,5]
        reference.genotypes.v0 <- cbind(input.map.v2, haplotypeToGenotype(reference.matrix, reference.genders))
    }


    # calculate allele frequencies form reference data
    if (is.null(reference.ped.map)) {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(input.genotypes.v0[,5:ncol(input.genotypes.v0)]), input.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,5:ncol(input.genotypes.v0)])
    } else {
        pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(reference.genotypes.v0[,5:ncol(reference.genotypes.v0)]), reference.ped.v1[,5])
        input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,c(5:ncol(input.genotypes.v0))])
    }
    colnames(input.genotypes.v1) <- c("chr", "snp_id", "pos_M","pos_bp", "freq", isolate.names)
    cat(paste("Begin filtering of ",length(isolate.names)," isolates and ",nrow(input.genotypes.v1)," SNPs...\n",sep=""))


    # remove SNPs with low population MAF
    input.genotypes.v2 <- subset(input.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf) #removing SNPs with AF>0.99 and AF<0.01
    if (nrow(input.genotypes.v2) == 0)
        stop("0 SNPs remain after MAF removal")
    cat(paste(nrow(input.genotypes.v2)," SNPs remain after MAF removal...\n",sep=""))


    # remove snps with high missingness
    snp.missingness    <- calculateMissingness(as.matrix(t(input.genotypes.v2[,6:ncol(input.genotypes.v2)])))
    input.genotypes.v3 <- input.genotypes.v2[snp.missingness <= snp.max.missing,]
    if (nrow(input.genotypes.v3) == 0)
        stop("0 SNPs remain after missingness removal")
    cat(paste(nrow(input.genotypes.v3)," SNPs remain after missingness removal...\n",sep=""))


    # remove samples with high missingness
    isolate.missingness <- round(calculateMissingness(as.matrix(input.genotypes.v3[,6:ncol(input.genotypes.v3)])),digits=3)
    if (length(isolate.names[isolate.missingness > isolate.max.missing]) > 0) {
        my.remove <- isolate.names[isolate.missingness > isolate.max.missing]
        warning("isolates removed due to genotype missingness: ",paste(my.remove, collapse=", "))
        sample.keep        <- input.ped.v1[isolate.missingness <= isolate.max.missing,1:6]
        input.genotypes.v4 <- input.genotypes.v3[,c(1:5, which(isolate.missingness <= isolate.max.missing) + 5)]
        if(nrow(sample.keep) < 1) stop(paste("All isolates removed with missingness > ",isolate.max.missing*100,"%. No isolates remaining.",sep=""))
    } else {
        sample.keep        <- input.ped.v1[,1:6]
        input.genotypes.v4 <- input.genotypes.v3
    }
    colnames(sample.keep) <- c("fid", "iid", "pid", "mid", "moi", "aff")
    if ((ncol(input.genotypes.v4)-5) == 0) {
        stop("0 samples remain after missingness removal")
    }
    cat(paste(ncol(input.genotypes.v4)-5," isolates remain after missingness removal...\n",sep=""))


    return.genotypes <- list(sample.keep, input.genotypes.v4)
    names(return.genotypes) <- c("pedigree", "genotypes")
    return(return.genotypes)

}


