## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
library(data.table)
test.data = "test.vcf.gz"
meta = "metadata.txt"
num.sample = as.numeric(system(sprintf("bcftools query -l %s | wc -l",test.data), intern = TRUE))
num.sites = as.numeric(system(sprintf("bcftools view -H %s | wc -l",test.data), intern = TRUE))
if(num.sample==0 || num.sites==0){
    stop("incorrect test vcf file")
}

metadata = fread(meta)
if(nrow(metadata)==0 || ncol(metadata)==0){
    stop("incorrect test metadata file")
}

