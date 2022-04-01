# rSNPdata     
R package to manage and analyse Plasmodium falciparum whole genome SNPs data.     

## requirements      
before installing this package, make sure the following tools and packages are installed:      
#### tools     
[bcftools](http://www.htslib.org/download/) version 1.9 or above   
[tabix](http://www.htslib.org/doc/tabix.html)  
[vcflib](https://github.com/vcflib/vcflib)

#### R packages     
`data.table`,`dplyr`,`SeqArray`,`moimix`,`statip`,`parallel`

## intallation   
``` {r}
library(devtools)
devtools::install_github("FellouMada/rSNPdata", build_vignettes = TRUE)
library(rSNPdata)
```

Short index:

- [Install](#Install)                                                             
- [DATA DESCRIPTION FUNCTIONS](#DATA DESCRIPTION)    
  * [get_snpdata](#DATA DESCRIPTION)
  * [print](#DATA DESCRIPTION)
  * [compute_MAF](#DATA DESCRIPTION)
  * [calculate_Fws](#DATA DESCRIPTION)
- [DATA FILTRATION FUNCTIONS](#DATA FILTRATION)
  * [filter_snps_samples](#DATA FILTRATION)
  * [select_chrom](#DATA FILTRATION)
  * [drop_snps](#DATA FILTRATION)
  * [drop_samples](#DATA FILTRATION)
- [DATA TRANSFORMATION FUNCTIONS](#DATA TRANSFORMATION)
  * [phase_mixed_genotypes](#DATA TRANSFORMATION)
  * [impute_missing_genotypes](#DATA TRANSFORMATION)

## Install   
``` {r eval=FALSE}
library(devtools)
devtools::install_github("FellouMada/rSNPdata", build_vignettes = TRUE)
library(rSNPdata)
```

## manual  
```{r}
browseVignettes("rSNPdata")
```

# DATA DESCRIPTION    

| function name | description |    
| :-------------- | :---------- |     
| [get_snpdata](./doc/get_snpdata.md) | Create SNPdata onject. the functions in this package require a `SNPdata` object. This can be generated with this function | 
| [print](./doc/print.md) | print the SNPdata object | 
| [compute_MAF](./doc/maf.md) | calculate the snp minor allele frequency based on the allelic depth and minor allele frequency | 
| [calculate_Fws](./doc/fws.md) | calculate the within-host genetic diversity index  | 

# DATA FILTRATION


# DATA TRANSFORMATION


