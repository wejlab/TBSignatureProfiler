## TBSignatureProfiler <img src="https://github.com/compbiomed/TBSignatureProfiler-docs/blob/master/lungs1.png?raw=true" align="right" width="170" />

[![](https://img.shields.io/github/last-commit/compbiomed/TBSignatureProfiler.svg)](https://github.com/compbiomed/TBSignatureProfiler/commits/master)
[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![codecov](https://codecov.io/gh/compbiomed/TBSignatureProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/TBSignatureProfiler)


## Documentation
Please check our website at https://compbiomed.github.io/TBSignatureProfiler-docs/.

## What is the TBSignatureProfiler?
The TBSignatureProfiler is an R package that allows researchers to easily profile RNA-seq data using signatures for TB progression, TB disease, and other TB disease states, using common signature profiling tools including ASSIGN, GSVA, and ssGSEA.

## Installation

The TBSignatureProfiler requires R 4.0 and Bioconductor 3.10.

* Install the development version of the package from Github 

``` r
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("compbiomed/TBSignatureProfiler")
```

* Install the release version of the package from Bioconductor

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TBSignatureProfiler")

```
