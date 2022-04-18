# NAME

**calculate_relatedness**

# USAGE

**calculate_relatedness**(snpdata, mat.name="GT", family="Location", sweepRegions=NULL, groups=c("Chogen","DongoroBa"))

# DESCRIPTION

**calculate_relatedness** is a function that measures the relatedmess between all pairs of isolates for any given pair of populations.

# OPTIONS
```
snpdata: an object of class SNPdata
mat.name: the name of the genotype table to be used. default="GT"
from: the name of the column, in the metadata table, to be used to represent the sample's population
sweepRegions: a data frame with the genomic coordinates of the regions of the genome to be discarded
groups: a vector of character. If specified, relatedness will be generated between these groups
```

# RETURN
a SNPdata object with an additional fied: relatedness. This is a list that contains the relatesness data frame and its correspondent matrix. The data frame will contain the following columns:
```
     1.  iid1: sample ID 1 in pair
     2.  iid2: sample ID 2 in pair
     3.  relatedness: the inferred relateness value
```
