# TBSignatureProfiler (development version)

## Bug Fixes
* Fixed some package tests for the profile.R, OriginalModel.R scripts to increase coverage

## Major Changes
* Added the gene_update parameter to runTBsigProfiler() to allow users to check signature/sample gene names for excel mogrified or outdated gene symbols
* Updated gene symbols for original models functions
* Added Chendi_HIV_2 signature (8/24)
* Removed Mendelsoh_RISK_11 as it was originally published in Darboe_RISK_11. Changed Darboe_RISK_11 to RISK11 in common objects. (Thomas Scriba) 

## Minor Changes
* Rewrote runTBsigProfiler() to be shorter and easier to maintain

# TBSignatureProfiler 1.5.0

* 68 signatures currently available

## Bug Fixes
* Fixed incorrect names of signatures in `Tbcommon` and `common_sigAnnotData` objects

## Major Changes
* Added Zimmer_RES_3, Gong_OD_4, Bloom_RES_268, and Bloom_RES_558
* Added Sivakumaran_11 and Mendelsoh_RISK_11 signatures
* Added Estevez_133, Estevez_259, LauxdaCosta_OD_3, and Maertzdorf_15 signatures to the package
* Added Chen_HIV_4, Gliddon_HIV_3, Gliddon_2_OD_4, Kulkarni_HIV_2, and Heycken_FAIL_22 signatures
* Added a COVIDsignatures object to the package that can be used to profile COVID-19 gene transcript signatures, thanks to collaborator Dylan Sheerin (WEHI)
* Added functions to evaluate some signatures using their original models from Johnson lab member Xutao Wang

## Minor Changes
* Added mention of COVIDsignatures object to main vignette and website
* Included the OG models tutorial on website (Xutao Wang)
* Added `addTBsignature()` to more easily facilitate updating signatures in package
* Added pROC option to obtain confidence intervals on AUC values as part of \code{tableAUC()}
* Added citation for newly published paper

# TBSignatureProfiler 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.