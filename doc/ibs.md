# NAME

**calculate_IBS**

# USAGE

**calculate_IBS**(snpdata, mat.name="GT")

# DESCRIPTION

**calculate_IBS** is a function that generates the dissimilarity matrix (1-IBS where IBS is the identity by state that represents the proportion of similar loci between any given pair of samples)  

# OPTIONS
```
snpdata: an object of class SNPdata
mat.name: the name of the genotype table to be used for IBS calculation
```

# RETURN
a SNPdata object with an additional field: IBS. This is an N*N matrix (where N is the number of samples)


