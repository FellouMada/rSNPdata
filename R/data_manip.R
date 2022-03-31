
require("data.table")
require("dplyr")
require("SeqArray")
require("moimix")
require("statip")
require("parallel")

#' Generate the SNPdata object
#' 
#' This function generate the input data needed for whole genome SNP data genotyped from malaria parasite 
#' @param vcf.file the input VCF file
#' @param meta.file the metadata file
#' @param output.dir the path to the folder where to store the output files
#' @return an object of class SNPdata with 5 elements
#' \describe{
#'   \item{meta}{data frame that contains the sample metadata}
#'   \item{details}{a data frame with SNPs genomic coordinates, fraction of missing data per SNP}
#'   \item{GT}{an integer matrix with the genotype data. 0='reference allele', 1='alternate allele', 2='mixed allele', NA='missing allele'}
#'   \item{vcf}{the vcf file from which the data is generated}
#'   \item{index}{an integer}
#'   }
#' @details use the print(snpdata) function to print the created object
#' @usage snpdata=get_snpdata(vcf.file = "file.vcf.gz", meta.file = "file.txt", output.dir = "/path/to/output/dir")
#' @export

get_snpdata = function(vcf.file=NA, meta.file=NA, output.dir=NA){
    if(is.na(vcf.file)) stop("Please provide an input VCF file")
    if(!file.exists(vcf.file)) stop(vcf.file, "not found!")
    if(is.na(meta.file)) stop("Please provide a metadata file")
    if(!file.exists(meta.file)) stop(meta.file, "not found!")
    if(is.na(output.dir)) stop("Please provide an output directory")
    if(!dir.exists(output.dir)) stop(output.dir, "not found!")
    
    ## get the sample IDs
    ids = paste0(output.dir,'/','SampleIDs.txt')
    system(sprintf("bcftools query -l %s > %s", vcf.file, ids))
    sampleIDs = fread(ids, header = FALSE)
    
    ## extracting the good quality SNPs 
    # filtered = paste0(output.dir,'/','Filtered.vcf.gz')
    # system(sprintf("bcftools view --threads 4 -i'N_ALT=1 && TYPE=\"%s\" && MQ>=%d && FORMAT/DP>=%d && FILTER=\"PASS\" && CDS && VQSLOD>=2' -o %s -Oz %s", variant, min.mq, min.dp, filtered, vcf.file))  #&& VQSLOD>=3
    
    ## extracting the genotype data
    genotypes = paste0(output.dir,'/','Genotypes.txt')
    expression = '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n'
    system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf.file, genotypes))
    genotypeF = fread(genotypes, header = FALSE, nThread = 4)
    names(genotypeF) = c("Chrom","Pos","Ref","Alt","Qual",sampleIDs$V1)
    
    ## the tables
    details = genotypeF %>% select(Chrom,Pos,Ref,Alt,Qual)
    names(sampleIDs) = "sample"
    snps = as.matrix(subset(genotypeF, select=-c(1:5)))
    snps[snps=="0/0"]="0"
    snps[snps=="1/1"]="1"
    snps[snps=="0/1" | snps=="1/0"]="2"
    snps[snps=="./." | snps==".|."]=NA
    snps=apply(snps, 2, function(x) as.integer(x))
    meta = add_metadata(sampleIDs,meta.file)
    meta$percentage.missing.sites = colSums(is.na(snps))/nrow(snps)
    details$percentage.missing.samples = rowSums(is.na(snps))/ncol(snps)
    
    snp.table = list(meta, details, snps, vcf.file, index=0)
    names(snp.table) = c("meta","details","GT","vcf","index")
    class(snp.table)="SNPdata"
    snp.table
}

add_metadata = function(sampleIDs, metadata){
    meta = fread(metadata, key = "sample")
    setkey(sampleIDs, "sample")
    samples=meta$sample
    if(any(!(samples %in% sampleIDs$sample))){
        warning("Incomplete meta data - should include all samples")
    }
    meta=data.frame(sample=samples) %>% left_join(meta,by="sample")
    meta
}

#' Print SNPdata
#' @param snpdata SNPdata object
#' @return NULL
#' @export
print.SNPdata=function(snpdata){
    print(head(snpdata$meta))
    print(head(snpdata$details))
    cat(sprintf("Data contains: %d samples for %d snp loci\n",dim(snpdata$GT)[2],dim(snpdata$GT)[1]))
    cat(sprintf("Data is generated from: %s \n",snpdata$vcf))
}

#' Filter loci and samples (requires bcftools and tabix to be installed)
#' 
#' This function generate the input data needed for whole genome SNP data genotyped from malaria parasite 
#' @param snpdata a SNPdata object
#' @param min.qual the minimum call quality score below which a loci will be discarded. default=10
#' @param max.missing.sites the maximum fraction of missing sites above which a sample should be discarded. default=0.2  
#' @param max.missing.samples the maximum fraction of missing samples above which a loci should be discarded. default=0.2  
#' @param maf.cutoff the MAF cut-off. loci with a MAF < maf.cutoff will be discarded
#' @return a SNPdata object
#' @usage snpdata = filter_snps_samples(snpdata, min.qual=10, max.missing.sites=0.2, max.missing.samples=0.2, maf.cutoff=0.01)
#' @export
filter_snps_samples=function (snpdata, min.qual=10, max.missing.sites=0.2, max.missing.samples=0.2, maf.cutoff=0.01){
    x=snpdata$details
    fields = c("GT","Phased","Phased_Imputed")
    if (missing(min.qual) & missing(max.missing.sites) & missing(max.missing.samples)) 
        return(snpdata)
    else {
        idx = which(x$Qual >= min.qual & x$percentage.missing.samples <= max.missing.samples & x$MAF >= maf.cutoff)
        if(length(idx)>0 & length(idx)<nrow(snpdata$details)){
            x = x[idx,]
            snpdata$details = x
            for(field in fields){
                if(field %in% names(snpdata)){
                    snpdata[[field]] = snpdata[[field]][idx,]
                }
            }
            f2c = x %>% select(Chrom, Pos)
            output.dir=dirname(snpdata$vcf)
            fwrite(f2c, paste0(output.dir,"/loci_to_be_retained.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = 4)
            snpdata$vcf = remove_snps_from_vcf(snpdata$vcf, "loci_to_be_retained.txt", output.dir, index = snpdata$index)
        }else if(length(idx)==0){
            stop("No locus in VCF file has satisfied specified the QC metrics")
        }else if(length(idx)==nrow(snpdata$details)){
            message("all loci have satisfied the specified QC metrics")
        }
        
        idx = which(snpdata$meta$percentage.missing.sites<=max.missing.sites)
        if(length(idx)>0 & length(idx)<nrow(snpdata$meta)){
            cat("the following samples will be removed: \n",paste(snpdata$meta$sample,collapse = "\n"))
            snpdata$meta = snpdata$meta[idx,]
            # x = x[,idx]
            fwrite(snpdata$meta$sample, paste0(output.dir,"/samples_to_be_dropped.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = 4)
            snpdata$vcf = remove_samples_from_vcf(snpdata$vcf, "samples_to_be_dropped.txt", output.dir, index = snpdata$index)
        }else if(length(idx)==0){
            stop("No sample in VCF file has satisfied the specified QC metrics")
        }else if(length(idx)==nrow(snpdata$meta)){
            message("all samples have satisfied the specified QC metrics")
        }
    }
    
    snpdata$index=snpdata$index+1
    snpdata
}

#' Calculate minor allele frequency (MAF)
#' 
#' Uses the SNPdata object to calculate the MAF at every loci 
#' @param snpdata a SNPdata object
#' @param include.het whether to account for the heterozygous allele or not. this is only used when mat.name="GT"
#' @param mat.name the name of the matrix to use. default is "GT"
#' @return a SNPdata object with 2 additional columns in the details data frame 
#' \describe{
#'   \item{MAF}{minor allele frequency of each snps}
#'   \item{MAF_allele}{1 if the alternate allele is the minor allele. 0 otherwise}
#' }
#' @details if include.het=FALSE, the mixed allele will not be considered in the MAF calculation
#' @usage snpdata = compute_MAF(snpdata, include.het=FALSE, mat.name="GT")
#' @export
compute_MAF = function(snpdata, include.het=TRUE, mat.name="GT"){
    x = snpdata[[mat.name]]
    ref = rowSums(x==0, na.rm = TRUE)
    alt = rowSums(x==1, na.rm = TRUE)
    het = rowSums(x==2, na.rm = TRUE)
    if(!include.het){
        res = apply(cbind(ref,alt), 1, getMaf)
    }else{
        res = apply(cbind(ref,alt,het), 1, getMaf)
    }
    if(!("MAF" %in% names(snpdata$details))){
        snpdata$details$MAF=as.numeric(res[1,])
        snpdata$details$MAF_allele = as.factor(as.character(as.numeric(round(res[2,]))))
        levels(snpdata$details$MAF_allele) = dplyr::recode_factor(snpdata$details$MAF_allele, REF="0", ALT="1", HET="2", REF_ALT="3",REF_ALT_HET="4")
    }else{
        new.maf = paste0("MAF_",mat.name)
        snpdata$details[[new.maf]] = as.numeric(res[1,])
    }
    snpdata
}

getMaf = function(mat){
    if(length(mat)==2){
        if(mat[1]<mat[2]){
            maf = mat[1]/sum(mat[1],mat[2])
            allele = 0
        }else if(mat[1]>mat[2]){
            maf = mat[2]/sum(mat[1],mat[2])
            allele = 1
        }else{
            maf = mat[2]/sum(mat[1],mat[2])
            allele = 3
        }
    }else{
        if(mat[1]<mat[2]){
            minor = mat[1]
            allele = 0
        }else if(mat[1]>mat[2]){
            minor = mat[2]
            allele = 1
        }else{
            minor = mat[2]
            allele = 3
        }
        
        if(minor<mat[3]){
            maf = minor/sum(mat[1],mat[2],mat[3])
            allele = 3
        }else if(minor>mat[3]){
            maf = mat[3]/sum(mat[1],mat[2],mat[3])
            allele = 2
        }else{
            maf = mat[3]/sum(mat[1],mat[2],mat[3])
            allele = 4
        }
    }
    return(c(maf, allele))
}


#' Calculate the complexity of infection in every sample (Fws)
#' 
#' Fws is the within host genetic diversity
#' @param snpdata a SNPdata object
#' @return a SNPdata object with an additional column in the meta table 
#' \describe{
#'   \item{Fws}{within host genetic diversity value}
#' }
#' @usage snpdata =  calculate_Fws(snpdata)
#' @export
calculate_Fws = function(snpdata){
    vcf = snpdata$vcf
    gdsFile = paste0(dirname(vcf),'/','data.gds')  #
    seqVCF2GDS(vcf, gdsFile)  
    my_vcf = seqOpen(gdsFile)     
    # seqSummary(my_vcf)
    sample.id = seqGetData(my_vcf, "sample.id")     
    coords = getCoordinates(my_vcf)      
    seqSetFilter(my_vcf, variant.id = coords$variant.id[coords$chromosome != "Pf3D7_API_v3"])     
    seqSetFilter(my_vcf, variant.id = coords$variant.id[coords$chromosome != "Pf_M76611"])      
    fws_overall = getFws(my_vcf)      #estimate the MOI
    fws_overall = data.table(cbind(as.character(sample.id), as.numeric(fws_overall)))
    names(fws_overall) = c("sample","Fws")
    meta=data.frame(snpdata$meta) %>% left_join(fws_overall,by="sample")
    snpdata$meta = meta
    snpdata
}

#' Phase mixed genotypes
#' 
#' mixed genotype phasing based on the number of reads supporting each allele of the heterozygous site. Simulation is performed 100 times 
#' @param snpdata a SNPdata object
#' @return a SNPdata object with an additional table named as "Phased". this will contain the phased genotypes 
#' @details when both alleles are not supported by any read or the total number of reads supporting both alleles at a given site is < 5, the genotype will be phased based on a bernouli distribition using the MAF as a parameter. Similarly, when the total number of reads is > 5 and the number of reads supporting one of the allele is not 2 times the number of the other allele, the genotype is phased using a bernouli distribution
#' @usage snpdata = phase_mixed_genotypes(snpdata)
#' @export
phase_mixed_genotypes = function(snpdata){
    vcf = snpdata$vcf
    expression = '%CHROM\t%POS[\t%AD]\n'
    tmp = paste0(dirname(vcf),"/tmp")
    system(sprintf("mkdir -p %s",tmp))
    ad = paste0(tmp,'/AllelicDepth.txt')
    system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, ad))
    depth = fread(ad, nThread = 4)
    depth = as.matrix(subset(depth, select = -c(1:2)))
    path = paste0(dirname(vcf),"/phasing")
    system(sprintf("mkdir -p %s", path))
    correlations = numeric(length = 100)
    for(i in 1:100){
        cat("running simulation ",i,"\n")
        tmp.snpdata = snpdata
        mat = apply(tmp.snpdata$GT, 1, phaseData, depth=depth)
        tmp.snpdata[["Phased"]]=t(mat)
        saveRDS(t(mat), paste0(path,"/sim",i,".RDS"))
        res.snpdata = compute_MAF(tmp.snpdata, include.het=FALSE, mat.name="Phased")
        correlations[i] = cor(res.snpdata$details[["MAF_Phased"]], res.snpdata$details[["MAF"]])
    }
    idx = which(correlations==max(correlations,na.rm = TRUE))
    snpdata[["Phased"]] = readRDS(paste0(path,"/sim",idx[1],".RDS"))
    system(sprintf("rm -rf %s", path))
    snpdata
}

phaseData = function(genotype, depth){
    idx = which(genotype==2)
    
    for(j in idx){
        ref = as.numeric(unlist(strsplit(depth[j],','))[1])
        alt = as.numeric(unlist(strsplit(depth[j],','))[2]) 
        if(ref==0 & alt==0){
            ref.count=sum(genotype==0, na.rm = TRUE)
            alt.count=sum(genotype==1, na.rm = TRUE)
            if(ref.count<alt.count) genotype[j] = 0
            else if(ref.count>alt.count) genotype[j] = 1
            else genotype[j] = rbern(1, ref.count/(ref.count+alt.count))
        }else if(ref!=0 & alt!=0){
            if((ref+alt)>=5 & (ref>=(2*alt) | alt>=(2*ref))){
                if(ref<alt) genotype[j] = 0
                else if(ref>alt) genotype[j] = 1
                else{
                    ref.count=sum(genotype==0, na.rm = TRUE)
                    alt.count=sum(genotype==1, na.rm = TRUE)
                    if(ref.count<alt.count) genotype[j] = 0
                    else if(ref.count>alt.count) genotype[j] = 1
                    else genotype[j] = rbern(1, ref.count/(ref.count+alt.count))
                }
            }else{
                ref.count=sum(genotype==0, na.rm = TRUE)
                alt.count=sum(genotype==1, na.rm = TRUE)
                if(ref.count<alt.count) genotype[j] = 0
                else if(ref.count>alt.count) genotype[j] = 1
                else genotype[j] = rbern(1, ref.count/(ref.count+alt.count))
            }
        }else if(ref==0 | alt==0){
            if(ref==0 & alt>=5) genotype[j] = 1
            else if(ref==0 & alt<5) genotype[j] = rbern(1, alt.count/(ref.count+alt.count))
            if(alt==0 & ref>=5) genotype[j] = 0
            else if(alt==0 & ref<5) genotype[j] = rbern(1, ref.count/(ref.count+alt.count))
        }
    }
    genotype
}


#' Impute missing genotypes
#' 
#' missing genotype imputation based on the MAF at any given locus 
#' @param snpdata a SNPdata object
#' @param genotype the genotype table from which the missing data will be imputed
#' @return a SNPdata object with an additional table named as "Phased_Imputed" 
#' @details when both alleles are not supported by any read or the total number of reads supporting both alleles at a given site is < 5, the genotype will be phased based on a bernouli distribition using the MAF as a parameter. Similarly, when the total number of reads is > 5 and the number of reads supporting one of the allele is not 2 times the number of the other allele, the genotype is phased using a bernouli distribution
#' @usage snpdata = impute_missing_genotypes(snpdata)
#' @export
impute_missing_genotypes = function(snpdata, genotype="GT"){
    cat("the mixed genotypes will be phased from ",genotype," table\n")
    field=genotype
    # genotype = snpdata[[genotype]]
    path = paste0(dirname(snpdata$vcf),"/imputing")
    system(sprintf("mkdir -p %s", path))
    correlations = numeric(length = 100)
    for(i in 1:100){
        cat("running simulation ",i,"\n")
        tmp.snpdata = snpdata
        mat = apply(tmp.snpdata[[field]], 1, impute)
        tmp.snpdata[["Imputed"]]=t(mat)
        saveRDS(t(mat), paste0(path,"/sim",i,".RDS"))
        res.snpdata = compute_MAF(tmp.snpdata, include.het=FALSE, mat.name="Imputed")
        correlations[i] = cor(res.snpdata$details[["MAF_Imputed"]], res.snpdata$details[["MAF"]])
    }
    idx = which(correlations==max(correlations,na.rm = TRUE))
    snpdata[["Imputed"]] = readRDS(paste0(path,"/sim",idx[1],".RDS"))
    system(sprintf("rm -rf %s", path))
    snpdata
}

impute = function(genotype){
    idx = which(genotype==2)
    for(j in idx){
        ref = length(which(genotype==0))
        alt = length(which(genotype==1))
        if(ref<alt) maf=ref/(ref+alt)
        else maf=alt/(ref+alt)
        genotype[j] = rbern(1, maf)
    }
    genotype
}

#' Select data from specified chromosomes
#' 
#' return data for specified chromosomes only 
#' @param snpdata a SNPdata object
#' @param chrom a vector of chromosomes
#' @return a SNPdata object with only the data from the specified chromosomes
#' @usage  chrom_snpdata = select_chrom(snpdata, chrom="Pf3D7_07_v3")
#' @export
select_chrom = function(snpdata, chrom="all"){
    system(sprintf("tabix %s", snpdata$vcf))
    m = which(names(snpdata) %in% c("meta","vcf","index"))
    fields = names(snpdata)[-m]
    if(chrom=="all"){
        return(snpdata)
    }
    res = list()
    for(chr in chrom){
        chrom.snpdata = snpdata
        idx = which(chrom.snpdata$details$Chrom==chr)
        for(field in fields){
            res[[field]] = rbind(res[[field]], chrom.snpdata[[field]][idx,])
        }
    }
    chrom.vcf = paste0(dirname(chrom.snpdata$vcf),"/target_chrom.vcf.gz")
    if(file.exists(chrom.vcf)){
        system(sprintf("rm -f %s", chrom.vcf))
    }
    if(length(chrom)>1){
        tmp.xme = paste0(dirname(chrom.snpdata$vcf),"/target_chrom.txt")
        fwrite(chrom, tmp.xme, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        system(sprintf("bcftools view -R %s %s -o %s -O z", tmp.xme, chrom.snpdata$vcf, chrom.vcf))
    }else{
        system(sprintf("bcftools view -r\"%s\" %s -o %s -O z", chrom, chrom.snpdata$vcf, chrom.vcf))
    }
    res$vcf = chrom.vcf
    res$meta = snpdata$meta
    res
}

#' Drop set a SNPs 
#' 
#' remove a set of SNPs from the SNPdata object
#' @param snpdata a SNPdata object
#' @param snp.to.be.dropped a data frame with 2 columns "Chrom" and "Pos"
#' @return a SNPdata object where the specified SNPs have been removed
#' @usage snpdata = drop_snps(snpdata, snp.to.be.dropped)
#' @export
drop_snps = function(snpdata, snp.to.be.dropped){
    if((is.data.frame(snp.to.be.dropped)) & (names(snp.to.be.dropped)%in%c("Chrom","Pos"))){
        idx = which(snpdata$details$Chrom%in%snp.to.be.dropped$Chrom & snpdata$details$Pos%in%snp.to.be.dropped$Pos)
        m = which(names(snpdata) %in% c("meta","vcf","index"))
        fields = names(snpdata)[-m]
        for(field in fields){
            snpdata[[field]] = snpdata[[field]][-idx]
        }
        f2c = snpdata$details %>% select(Chrom,Pos)
        tmp.file = paste0(dirname(snpdata$vcf),"/tmp.txt")
        fwrite(f2c, tmp.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", nThread = 4)
        snpdata$vcf = remove_snps_from_vcf(snpdata$vcf, "tmp.txt", path=dirname(snpdata$vcf), index=snpdata$index)
    }else{
        stop("value for 'snp.to.be.dropped' argument should a data frame whith Chrom and Pos columns.")
    }
    snpdata$index = snpdata$index+1
    snpdata
}

remove_snps_from_vcf = function(vcf, loci_to_be_retained, path, index=1){
    target.loci = paste0(path,"/",loci_to_be_retained)
    header = paste0(path,'/','Header.txt')
    body = paste0(path,'/','Body.txt')
    correctRows = paste0(path,'/','Good_snps.txt')
    filteredVcf = paste0(path,'/','Filtered_snps_',index,'.vcf')
    
    system(sprintf("bcftools view -h %s > %s", vcf, header))
    system(sprintf("bcftools view -H %s > %s", vcf, body))
    system(sprintf("awk -F'\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' %s %s > %s",target.loci,body,correctRows))
    system(sprintf("cat %s %s > %s", header, correctRows, filteredVcf))
    system(sprintf("bgzip %s", filteredVcf))
    system(sprintf("rm -f %s %s %s", header, body, correctRows))
    return(as.character(filteredVcf))
}

#' Drop samples  
#' 
#' remove a set of samples from the SNPdata object
#' @param snpdata a SNPdata object
#' @param samples.to.be.dropped a vector of samples to be dropped
#' @return a SNPdata object where the specified samples have been removed
#' @usage snpdata = drop_samples(snpdata, samples.to.be.dropped)
#' @export
drop_samples = function(snpdata, samples.to.be.dropped){
    if(length(samples.to.be.dropped)==0 | (any(!(samples.to.be.dropped %in%snpdata$meta$sample)))){
        stop("no provided samples or provided samples not found!")
    }
    idx = match(samples.to.be.dropped, snpdata$meta$sample)
    tmp_meta= snpdata$meta
    tmp_meta = tmp_meta[-(idx),]
    snpdata$meta = tmp_meta
    m = which(names(snpdata) %in% c("details","vcf","meta","index"))
    fields = names(snpdata)[-m]
    for(field in fields){
        idx = match(samples.to.be.dropped, colnames(snpdata[[field]]))
        # print(idx)
        tmp_meta=snpdata[[field]]
        tmp_meta = tmp_meta[,-(idx)]
        snpdata[[field]] = tmp_meta
    }
    tmp.file = paste0(dirname(snpdata$vcf),"/tmp.txt")
    write.table(snpdata$meta$sample, tmp.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    snpdata$vcf = remove_samples_from_vcf(snpdata$vcf, "tmp.txt", path=dirname(snpdata$vcf), index=snpdata$index)
    snpdata$index = snpdata$index+1
    snpdata
}

remove_samples_from_vcf = function(vcf, samples.to.be.retained, path, index=1){
    target.samples = paste0(path,"/",samples.to.be.retained)
    post.qc = paste0(path,'/','Post_QC_',index,'.vcf.gz')
    system(sprintf("bcftools view -S %s %s -o %s -O z", target.samples, vcf, post.qc))
    system(sprintf("rm -f %s", vcf))
    system(sprintf("tabix %s", post.qc))
    return(as.character(post.qc))
}







