

#' Calculate population level allele frequency and the within-sample allele frequency
#'
#' the population level allele frequency and the within-sample allele frequency are required for the deconvolution process
#' @param snpdata a SNPdata object
#' @return a list of 2 elements:
#' \enumerate{
#'   \item PLAF: a data frame with 5 column: CHROM, POS, PLAF (population level allele frequency), NRAF (non-reference allele frequency), MAF(the minor allele frequency from the non-reference allele frequency)
#'   \item WSAF: a matrix with the within sample allele frequency
#' }
#' @usage pop.level.allele.frequency = populationLevelAlleleFrequency(snpdata)
#' @export
populationLevelAlleleFrequency = function(snpdata)
{
    vcf = snpdata$vcf
    expression = '%CHROM\t%POS[\t%AD]\n'
    tmp = paste0(dirname(vcf),"/tmp")
    system(sprintf("mkdir -p %s",tmp))
    ad = paste0(tmp,'/AllelicDepth.txt')
    system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, ad))
    depth = fread(ad, nThread = 4)
    f2c = depth %>% select(V1,V2)
    depth = as.matrix(depth%>% select(-c(V1,V2)))
    wsaf = matrix(-9, nrow = dim(depth)[1], ncol = dim(depth)[2])
    colnames(wsaf) = snpdata$meta$sample
    plaf = vector('numeric',dim(depth)[1])
    pb = txtProgressBar(min = 0, max = dim(wsaf)[1], initial = 0,style = 3, char = "*")
    for(i in 1:dim(wsaf)[1]){
        sumREF = 0; sumALT = 0
        for(j in 1:dim(wsaf)[2]){
            if(is.na(depth[i,j]) | depth[i,j]=="."){
                ref=alt=0
            }else{
                ref=as.integer(unlist(strsplit(depth[i,j],','))[1])
                alt=as.integer(unlist(strsplit(depth[i,j],','))[2])
            }
            sumREF=sumREF+ref
            sumALT=sumALT+alt
            if(alt==0 & ref==0)
                wsaf[i,j]=0
            else
                wsaf[i,j]=alt/(ref+alt)
        }
        plaf[i]=sumALT/(sumREF+sumALT)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    plaf = data.frame(cbind(f2c,plaf), stringsAsFactors = FALSE)
    names(plaf) = c('CHROM','POS','PLAF')
    non.ref.allele.frequency = rowMeans(wsaf, na.rm = TRUE, dims = 1)
    maf = vector("numeric",length(non.ref.allele.frequency))
    for(i in 1:length(maf))
        maf[i]=min(non.ref.allele.frequency[i], 1-non.ref.allele.frequency[i])
    plaf$NRAF = non.ref.allele.frequency
    plaf$MAF = maf
    system(sprintf("rm -rf %s", tmp))
    list(PLAF = plaf, WSAF=wsaf)
}



