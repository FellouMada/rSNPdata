# NAME

**filter_snps_samples**

# USAGE

**filter_snps_samples**(snpdata, min.qual=10, max.missing.sites=0.2, max.missing.samples=0.2, maf.cutoff=0.01)

# DESCRIPTION

**filter_snps_samples** is a function that is used to filter the SNPdata object. Note that, this must be used only after the MAF and Fws have been calculated

# OPTIONS
```
snpdata: an object of class SNPdata
min.qual: the minimum call quality of a SNP
max.missing.sites: the maximum fraction of missing sites. Any sample with a fraction of missing data > max.missing.sites will be discarded
max.missing.samples: the maximum fraction of missing samples. Any locus with a fraction of missing data > max.missing.samples will be discarded
maf.cutoff: the MAF cut-off. loci with MAF < maf.cutoff will be removed
```

# RETURN
a SNPdata object where poor quality samples and SNPs have been removed
