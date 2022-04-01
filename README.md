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
- [FUNCTIONS](#FUNCTIONS)    
  * [SNPdata object](#get SNPdata object)
  * [Metrics](#metrics)
  * [Phenotype](#phenotype)
  * [Genotype](#genotype)
  * [Transformation](#transformation)
  * [Statistics](#statistics)
  * [Scripts](#scripts)
- [Link library](#link-library)
- [Build from source](#build-from-source)
- [Development](#Development)
- [LICENSE](#LICENSE)
- [Credit work](#Credit)

## Install   
``` {r eval=FALSE}
library(devtools)
devtools::install_github("FellouMada/rSNPdata", build_vignettes = TRUE)
library(rSNPdata)
```

# FUNCTIONS     

## get SNPdata object    
| function name | description |    
| :-------------- | :---------- |     
| [get_snpdata](./doc/get_snpdata.md) | Create SNPdata onject. the functions in this package require a `SNPdata` object. This can be generated with this function| 


## manual  
```{r}
browseVignettes("rSNPdata")
```

