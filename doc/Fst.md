# NAME

**calculate_wcFst**

# USAGE

**calculate_wcFst**(snpdata, from="Country", groups=c("Senegal","Gambia"))

# DESCRIPTION

**calculate_wcFst** is a function that calculates the Weir & Cockerham's Fst. Note that the wcFst is calculated using the [vcflib[(https://github.com/vcflib/vcflib) tools

# OPTIONS
```
snpdata: an object of class SNPdata
from: the name of the column, in the metadata table, that contains the information about the groups the samples belong to
groups: a vector of characters. the elements should be taken from the `from` column
```

# RETURN
a SNPdata object with an additional fied, Fst. This is a list of data frames that contain the result from a specific comparison. Every data frame contains N rows (where N is the number of loci) and the following 7 columns:
```
     1.  Chrom: the chromosome ID
     2.  Pos: the SNPs position
     3.  AF_in_target: allele frequency in the first group
     4.  AF_in_background: allele frequency in the second group
     5.  wcFst: the resulting Fst values
     6.  wcFst_pvalue: the p-values associated to the Fst results
     7.  wcFst_Adj_pvalue_BH: the adjusted p-value for multiple testing using the Benjamini-Hochberg method
```

