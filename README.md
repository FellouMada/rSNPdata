# rSNPdata     
R package to manage and analyse Plasmodium falciparum whole genome SNPs data.     

## requirements      
before installing this package, make sure the following tools and packages are installed:      
#### tools     
[bcftools](http://www.htslib.org/download/)    
[tabix](http://www.htslib.org/doc/tabix.html)  
[vcflib](https://github.com/vcflib/vcflib)

#### R packages     
`data.table`,`dplyr`,`SeqArray`,`moimix`,`statip`,`parallel`

---
Index:

- [Install](#Install) 
- [Manual](#Manual)
- [DATA DESCRIPTION](#DESCRIPTION)    
  * [get_snpdata](#DESCRIPTION)
  * [print](#DESCRIPTION)
  * [compute_MAF](#DESCRIPTION)
  * [calculate_Fws](#DESCRIPTION)
- [DATA FILTRATION](#FILTRATION)
  * [filter_snps_samples](#FILTRATION)
  * [select_chrom](#FILTRATION)
  * [drop_snps](#FILTRATION)
  * [drop_samples](#FILTRATION)
- [DATA TRANSFORMATION](#TRANSFORMATION)
  * [phase_mixed_genotypes](#TRANSFORMATION)
  * [impute_missing_genotypes](#TRANSFORMATION)
---

## Install   
``` {r eval=FALSE}
library(devtools)
devtools::install_github("FellouMada/rSNPdata", build_vignettes = TRUE)
library(rSNPdata)
```

## Manual  
```{r}
browseVignettes("rSNPdata")
```

# DESCRIPTION    

| function name | description |    
| :-------------- | :---------- |     
| [get_snpdata](./doc/get_snpdata.md) | Create SNPdata onject. the functions in this package require a `SNPdata` object. This can be generated with this function | 
| [print](./doc/print.md) | print the SNPdata object | 
| [compute_MAF](./doc/maf.md) | calculate the snp minor allele frequency based on the allelic depth and minor allele frequency | 
| [calculate_Fws](./doc/fws.md) | calculate the within-host genetic diversity index  | 

# FILTRATION

| function name | description |    
| :-------------- | :---------- |     
| [filter_snps_samples](./doc/get_snpdata.md) | Create SNPdata onject. the functions in this package require a `SNPdata` object. This can be generated with this function | 
| [select_chrom](./doc/print.md) | print the SNPdata object | 
| [drop_snps](./doc/maf.md) | calculate the snp minor allele frequency based on the allelic depth and minor allele frequency | 
| [drop_samples](./doc/fws.md) | calculate the within-host genetic diversity index  | 


# TRANSFORMATION

| function name | description |    
| :-------------- | :---------- |     
| [phase_mixed_genotypes](./doc/get_snpdata.md) | Create SNPdata onject. the functions in this package require a `SNPdata` object. This can be generated with this function | 
| [impute_missing_genotypes](./doc/print.md) | print the SNPdata object | 


