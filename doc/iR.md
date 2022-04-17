# NAME

**calculate_iR**

# USAGE

**calculate_iR**(snpdata, mat.name="Phased", family="Location", number.cores=4)

# DESCRIPTION

**calculate_iR** is a function that calculates the iR index to determine loci with excess of IBD sharing. IBD here is calculated based on the [isoRelate](https://github.com/bahlolab/isoRelate) R package.

# OPTIONS
```
snpdata: an object of class SNPdata
mat.name: the name of the genotype table to be used. default is "Phased"
family: the name of the column, in the metadata table, to be used to represent the sample's population
number.cores: the number of cores to be used. default=4
```

# RETURN
a SNPdata object with an additional fied, iR. This is a list of data frames that contain the iR results from the specified pairs of populations. Every data frame contains N rows (where N is the number of loci) and the following 9 columns:
```
     1.  Chrom: the chromosome ID
     2.  snp_id: the SNPs identifiers
     3.  pos_M: the genetic map distance (centi morgans, cM)
     4.  pos_bp: the SNP position
     5.  pop: the population
     6.  subpop: the subpopulation
     7.  iR: the iR statistic
     8.  pvalue: the p-value associated to iR
     9.  adj_pvalue_BH: the adjusted p-value after multiple testing correction using the Benjamini-HochBerg method
     
```
