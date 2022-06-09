# rSNPdata     
R package to manage and analyse Plasmodium falciparum whole genome SNPs data.     

## requirements      
before installing this package, make sure the following tools are installed:      
#### tools     
[bcftools](http://www.htslib.org/download/)    
[tabix](http://www.htslib.org/doc/tabix.html)  
[vcflib](https://github.com/vcflib/vcflib)    
[vcftools](http://vcftools.sourceforge.net/)    
[BEDOPS](https://bedops.readthedocs.io/en/latest/index.html)

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
- [POPULATION GENETICS PARAMETERS](#PARAMETERS)
  * [calculate_wcFst](#PARAMETERS)
  * [calculate_LD](#PARAMETERS)
  * [calculate_IBS](#PARAMETERS)
  * [calculate_iR](#PARAMETERS)
  * [calculate_relatedness](#PARAMETERS)
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
| [filter_snps_samples](./doc/filter.md) | filter loci and samples from the SNPdata object | 
| [select_chrom](./doc/select_chrom.md) | select data for a provided list of chromosomes | 
| [drop_snps](./doc/remove_snp.md) | remove a set of SNPs from the SNPdata object | 
| [drop_samples](./doc/remove_sample.md) | remove a set of samples from the SNPdata object  | 


# TRANSFORMATION

| function name | description |    
| :-------------- | :---------- |     
| [phase_mixed_genotypes](./doc/phase.md) | Phase mixed genotypes based on number of read supporting each allele and Bernoulli distribution. This process will be repeated `nsim` times and data from iteration with highest correlation between MAF phased data and MAF raw data will be retained | 
| [impute_missing_genotypes](./doc/impute.md) | impute missing genotypes based on minor allele frequency and Bernoulli distribution. This process will be repeated `nsim` times and data from iteration with highest correlation between MAF imputed data and MAF raw data will be retained| 


# PARAMETERS

| function name | description |    
| :-------------- | :---------- |     
| [calculate_wcFst](./doc/Fst.md) | Calculate Weir & Cockerham's Fst. This is achieved using the vcflib tools | 
| [calculate_LD](./doc/ld.md) | Calculate LD R^2 between all pair of loci using vcftools| 
| [calculate_IBS](./doc/ibs.md) | Calculate identity by state matrix | 
| [calculate_iR](./doc/iR.md) | Calculate the iR index between pairs of populations to determine loci with excess of IBD sharing. The calculation is based on the isoRelate R package | 
| [calculate_relatedness](./doc/relatedness.md) | Calculate relatedness between every pair od isolates for all pairs of population. This is based on the model developed by [Aimee R. Taylor and co-authors](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009101) | 


