# NAME

**calculate_Fws**

# USAGE

**calculate_Fws**(snpdata)

# DESCRIPTION

**calculate_Fws** is a function that calculates the within-host genetic diversity (Fws). Note that the Fws is obtained using the [moimix](https://github.com/bahlolab/moimix) R package

# OPTIONS
```
snpdata: an object of class SNPdata
```

# RETURN
a SNPdata object with an additional column in the `meta` data frame:
```
     1.  Fws: the Fws across all samples
```
