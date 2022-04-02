# NAME

**compute_MAF**

# USAGE

**compute_MAF**(snpdata, include.het=FALSE, mat.name="GT")

# DESCRIPTION

**compute_MAF** is a function that calculates the minor allele frequency (MAF) for every locus from the specified genotype table specified in `mat.name` argument ("GT": calculate MAF from the raw data; "Phased": calculate MAF from the phased data; "Imputed": calculate MAF from the imputed data). When `include.het=TRUE`, the mixed genotypes will be considered in the MAF calculation

#OPTIONS
```
snpdata: an object of class SNPdata
include.het: whether to account for mixed genotypes or not
mat.name: the genotype matrix from which MAF will be calculated
```

# RETURN
a SNPdata object with 2 additional columns in the `details` data frame:
```
     1.  MAF: the minor allele frequency
     2.  MAF_allele: the least frequent allele  
```
