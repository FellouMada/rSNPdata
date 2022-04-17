# NAME

**calculate_LD**

# USAGE

**calculate_LD**(snpdata, min.r2=0.2, inter.chrom=FALSE, chroms=NULL)

# DESCRIPTION

**calculate_LD** is a function that is used to calculate LD R^2 between pairs of loci using the [vcftools](http://vcftools.sourceforge.net/) program. 

# OPTIONS
```
snpdata: an object of class SNPdata
min.r2: the minimum R^2 value below which the LD value is not reported
inter.chrom: whether to calculate inter-chromosomal LD
chroms: a vector of chromosome names for which LD should be calculated
```

# RETURN
a SNPdata object with an additional field: LD. this is a data frame with 6 columns:
```
     1.  CHR1: the chromosome of the first SNP
     2.  POS1: the position of the first SNP in the chromosome
     3.  CHR2: the chromosome of the second SNP
     4.  POS2: the position of the second SNP in the chromosome
     5.  N_INDV: the number of samples used to calculate the R^2
     6.  R^2: the resulting R^2 values
```
