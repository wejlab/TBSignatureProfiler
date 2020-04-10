## TBSignatureProfiler <img src="https://github.com/compbiomed/TBSignatureProfiler-docs/blob/master/lungs1.png?raw=true" align="right" width="170" />

[![Travis build status](https://travis-ci.org/compbiomed/TBSignatureProfiler.svg?branch=master)](https://travis-ci.org/compbiomed/TBSignatureProfiler)
[![codecov](https://codecov.io/gh/compbiomed/TBSignatureProfiler/branch/master/graph/badge.svg)](https://codecov.io/gh/compbiomed/TBSignatureProfiler)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## What is the TBSignatureProfiler?
The TBSignatureProfiler is an R package that allows researchers to easily profile RNA-seq data using signatures for TB progression, TB disease, and other TB disease states, using common signature profiling tools including ASSIGN, GSVA, and ssGSEA.

Tuberculosis (TB) is the leading cause of infectious disease mortality worldwide, causing, on average, nearly 1.4 million deaths per year. One of the most consistent issues faced in controlling TB stems from difficulty diagnosing the disease. Particular types of TB infections, including pediatric TB, cannot be accurately identified with traditional bacteria tests (e.g., via phlegm/sputum). As an alternative diagnosis mechanism, researchers have identified gene expression signatures in blood as effective disease biomarkers. A gene signature is defined as a set of genes with a characteristic pattern of gene expression which occurs as a result of a biological, cellular or molecular phenotype. To date, more than 30 TB signatures have been published by researchers, though most have relatively low cross-condition validation (e.g., testing TB in samples from diverse geographic and comorbidity backgrounds, where 2 or more chronic diseases are present). Furthermore, these signatures have never been formally collected and made available as a single unified resource.

The goal of the TBSignatureProfiler is to aggregate these signatures and create an efficient platform for their visual and quantitative comparison via an R package that characterizes gene signatures’ diagnostic ability in multiple comorbidity settings. The package uses computational methods to analyze expression levels in a sample for a given set of genes (that is, a gene expression profile). These profiles may then be scored to quantify gene expression patterns and correlated to TB outcomes, including distinguishing between active TB disease and latent TB infection (LTBI). In addition to profiling, the software enables gene signature accuracy prediction by bootstrapping area under the curve estimates (AUC, a statistical evaluation metric), and logistic regression with leave-one-out cross-validation (LOOCV). Visualization tools include sample-signature score heatmaps, bootstrap AUC and LOOCV boxplots, and results tables. An R shiny app is included that presents a graphical interface with which to use the package.

## Installation

TBSignatureProfiler is under development. You can install the devel version via
GitHub:

``` r
# install.packages("devtools")
devtools::install_github("compbiomed/TBSignatureProfiler")
```

## Run TBSignatureProfiler shiny app

``` r
library(TBSignatureProfiler)
TBSPapp()
```