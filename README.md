## TBSignatureProfiler <img src="https://github.com/wejlab/TBSignatureProfiler-docs/blob/master/lungs1.png?raw=true" align="right" width="170" />

[![](https://img.shields.io/badge/bioconductor-3.18-blue)](http://www.bioconductor.org/packages/release/bioc/html/TBSignatureProfiler.html)
[![codecov](https://codecov.io/gh/wejlab/TBSignatureProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/wejlab/TBSignatureProfiler)
[![build](http://www.bioconductor.org/shields/build/release/bioc/TBSignatureProfiler.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/TBSignatureProfiler/)

## What is the TBSignatureProfiler?
The TBSignatureProfiler (TBSP) is an easy-to-use R package for analysis of RNA-seq data using gene signatures for tuberculosis disease presence, risk, progression, treatment failure, and other states. In-package signature profiling is available using common gene set enrichment tools that include GSVA, singscore, and ssGSEA.

Users can analyze RNA-seq data across 70+ published gene signatures to compare signature performance with clear visualizations (e.g., heatmaps and ROC plots), specificity/sensitivity estimates for predicting binary variables, and cross-validated logistic regression.

The TBSP can be used as a standalone software package. Data exploration is also available with the accompanying interactive R Shiny application. The R Shiny app is currently functional but more limited in its capabilities than the command line form of the package.

## Documentation

Please check our website at https://wejlab.github.io/TBSignatureProfiler-docs/.

## Citation

We used the TBSP to compare 45 gene signatures and differentiate active TB from LTBI in malnourished individuals from India. You can read our paper here: [“Comparing tuberculosis gene signatures in malnourished individuals using the TBSignatureProfiler”](
https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-020-05598-z)

Please cite:

**Johnson, W.E., Odom, A., Cintron, C. et al. Comparing tuberculosis gene signatures in malnourished individuals using the TBSignatureProfiler. BMC Infect Dis 21, 106 (2021).**

## Installation

The TBSignatureProfiler requires R Version 4.1.

* Install the development version of the package from Github:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("wejlab/TBSignatureProfiler")
```

* Install the release version of the package from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TBSignatureProfiler")

```

* Run the TBSP R Shiny app!
``` r
library(TBSignatureProfiler)
TBSPapp()
```

