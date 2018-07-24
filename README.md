enrichedRanges [![Build Status](https://travis-ci.org/leekgroup/enrichedRanges.svg?branch=master)](https://travis-ci.org/leekgroup/enrichedRanges)
==============

Identify enrichment between two sets of genomic ranges.

# Installation instructions

Get R 3.1 or newer from [CRAN](http://cran.r-project.org/).

```S
## If needed
install.packages("devtools")

## Pre-requisites from Bioconductor
install.packages("BiocManager")
BiocManager::install(c("biovizBase", "GenomicRanges", "seqbias"))
    
## enrichedRanges
library("devtools")
install_github("leekgroup/enrichedRanges")
```


# Citation

Below is the citation output from using `citation("enrichedRanges")` in R. 
Please run this yourself to check for any updates on how to cite 
__enrichedRanges__.

---

To cite package __enrichedRanges__ in publications use:

Leonardo Collado-Torres and Andrew Jaffe (2014). enrichedRanges: Identify enrichment between two sets of genomic ranges. R package version 0.0.2.

A BibTeX entry for LaTeX users is

```
@Manual{,
    title = {enrichedRanges: Identify enrichment between two sets of genomic ranges},
    author = {Leonardo Collado-Torres and Andrew Jaffe},
    note = {R package version 0.0.2},
}
```

## Travis CI

This package is automatically tested thanks to [Travis CI](travis-ci.org) and [r-travis](https://github.com/craigcitro/r-travis). If you want to add this to your own package use:

```R
## Use devtools to create the .travis.yml file
library('devtools')
use_travis('yourPackage')

## Read https://github.com/craigcitro/r-travis/wiki to configure .travis.yml appropriately

## Add a status image by following the info at http://docs.travis-ci.com/user/status-images/
```
