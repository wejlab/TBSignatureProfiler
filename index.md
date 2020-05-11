<img src="https://github.com/compbiomed/TBSignatureProfiler-docs/blob/master/lungs1.png?raw=true" align="right" width="170" />

[![](https://img.shields.io/badge/bioconductor-3.11-3a6378.svg)](http://www.bioconductor.org/packages/release/bioc/html/TBSignatureProfiler.html)
[![codecov](https://codecov.io/gh/compbiomed/TBSignatureProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/TBSignatureProfiler)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue)](https://www.tidyverse.org/lifecycle/#maturing)
[![build](http://www.bioconductor.org/shields/build/release/bioc/TBSignatureProfiler.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/TBSignatureProfiler/)


## What is the TBSignatureProfiler?
The TBSignatureProfiler is an R package that allows researchers to easily profile RNA-seq data using signatures for TB progression, TB disease, and other TB disease states, using common signature profiling tools including ASSIGN, GSVA, and ssGSEA.

Tuberculosis (TB) is the leading cause of infectious disease mortality worldwide, causing, on average, nearly 1.4 million deaths per year. One of the most consistent issues faced in controlling TB stems from difficulty diagnosing the disease. Particular types of TB infections, including pediatric TB, cannot be accurately identified with traditional bacteria tests (e.g., via phlegm/sputum). As an alternative diagnosis mechanism, researchers have identified gene expression signatures in blood as effective disease biomarkers. A gene signature is defined as a set of genes with a characteristic pattern of gene expression which occurs as a result of a biological, cellular or molecular phenotype. To date, more than 30 TB signatures have been published by researchers, though most have relatively low cross-condition validation (e.g., testing TB in samples from diverse geographic and comorbidity backgrounds, where 2 or more chronic diseases are present). Furthermore, these signatures have never been formally collected and made available as a single unified resource.

The goal of the TBSignatureProfiler is to aggregate these signatures and create an efficient platform for their visual and quantitative comparison via an R package that characterizes gene signatures’ diagnostic ability in multiple comorbidity settings. The package uses computational methods to analyze expression levels in a sample for a given set of genes (that is, a gene expression profile). These profiles may then be scored to quantify gene expression patterns and correlated to TB outcomes, including distinguishing between active TB disease and latent TB infection (LTBI). In addition to profiling, the software enables gene signature accuracy prediction by bootstrapping area under the curve estimates (AUC, a statistical evaluation metric), and logistic regression with leave-one-out cross-validation (LOOCV). Visualization tools include sample-signature score heatmaps, bootstrap AUC and LOOCV boxplots, and results tables. An R shiny app is included that presents a graphical interface with which to use the package.

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
## Run TBSignatureProfiler shiny app

``` r
library(TBSignatureProfiler)
TBSPapp()
```