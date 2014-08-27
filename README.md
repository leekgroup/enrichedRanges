enrichedRanges
==============

Identify enrichment between two sets of genomic ranges.

# Installation instructions

Get R 3.1 or newer from [CRAN](http://cran.r-project.org/).

```S
## If needed
install.packages("devtools")

## Pre-requisites from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("biovizBase", "GenomicRanges", "seqbias"))
    
## enrichedRanges
library("devtools")
install_github("lcolladotor/enrichedRanges")
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
