# NAME

**impute_missing_genotypes**

# USAGE

**impute_missing_genotypes**(snpdata, genotype="Phased")

# DESCRIPTION

**impute_missing_genotypes** is a function that imputes the missing genotypes. Missing data can be phased from raw data (genotype="GT") or phased data (genotype="Phased")

# OPTIONS
```
snpdata: an object of class SNPdata
genotype: the name of the genotype matrix to be used
```

# RETURN
a SNPdata object with an additional field
```
     Imputed: the matrix of imputed data
```
