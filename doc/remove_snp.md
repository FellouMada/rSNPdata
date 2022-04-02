# NAME

**drop_snps**

# USAGE

**drop_snps**(snpdata, snp.to.be.dropped=NA, chrom=NA, start=NA, end=NA)

# DESCRIPTION

**drop_snps** is a function that removes a set of snps from the SNPdata object. the genomic coordinates of the loci to be removed could be given as a data frame with 2 columns: "Chrom" and "Pos" (in this case, the rest of the argument can be ignored or set to NA) or as a genomic region: chrom="Pf3D7_10_v3", start=100, end=500 (in this case, the snp.to.be.dropped argument can be ignored or set to NA)  

#OPTIONS
```
snpdata: an object of class SNPdata
snp.to.be.dropped: a data frame with 2 columns: Chrom and Pos
chrom: the chromosome to which to loci to be removed belong 
start: the starting position of the loci to be removed
end: the end position of the loci to be removed
```

# RETURN
a SNPdata object where the corresponding loci of the specified genomic coordinates have been removed
