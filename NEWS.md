# TBSignatureProfiler 1.19.0

## Minor Changes
* Updated gene name update mechanism in `runTBsigProfiler` to sum genes with identical names

## Bug Fixes
* Fixed issue with GSVA needing a non-data.frame object in `runTBsigProfiler`

# TBSignatureProfiler 1.17.0

## Minor Changes
* Added Li_3 signature
* Adjusted function to create documentation in Signature Addition vignette
* Fixed Zhao_Nano_6 signature documentation

## Bug Fixes
* Changed out .data pronouns to use proper tidyr variable calling

# TBSignatureProfiler 1.16.0

## Bug Fixes
* Removed references to compbiomed repository and replaced with wejlab.
* Fixed calls to GSVA package dependency to reflect new implementation (runTBsigProfiler)
* Updated vignette to reflect more recent results (using more signatures than when package was first released)

# TBSignatureProfiler 1.14.0

## Bug Fixes
* Fixed bug for `signatureBoxplot` by removing quotations around variables. This issue was introduced with the newest ggplot2 version. Thanks to Arthur VanValkenburg for identifying the issue and solution.
* Removed Zimmer_RES_3 which was an identical signature as the previously published Sweeney_OD_3.
* Fixed error in example for `compare_algs` caused by the signature being used.

## Minor Changes
* Added Vargas_18 and Vargas_45 signatures (doi: 10.1371/journal.pcbi.1010770)

# TBSignatureProfiler 1.12.0

## Bug Fixes
* Removed pROC argument from `TableAUC()` as it was causing a bug where the upper CI was always
the same as the AUC point estimate, even with `pROC=FALSE` as the default setting.
* Changed source package for `rowSums()` and `rowMeans()` to `base` from `BiocGenerics` due to an error on the Bioconductor testing servers.

## Minor Changes
* Added Chen_5 signature

# TBSignatureProfiler 1.10.0

## Bug Fixes
* Fixed the tableAUC bootstrapped confidence interval to be the 2.5 and 97.5 percentiles instead of the 5 and 95 percentiles
* Fixed the upper CI value for the pROC/DeLong AUC CI method in the bootstrapAUC function

## Major Changes
* Changed tableAUC confidence interval default to bootstrapped instead of DeLong (pROC argument)

## Minor Changes
* Updated the github and website introductions.
* Added Natarajan_7, Kaul_3 signatures
* Added Francisco_OD_2, Kwan_186 signatures
* Shortened some example run times

# TBSignatureProfiler 1.8.0

## Bug Fixes
* Fixed gene in RESPONSE5 (PNN to RP11-295G20.2) in TBsignatures and TBcommon objects. (Stanley M. Kimbung)
* Fixed a bug in the singscore algorithm called from runTBsigProfiler() that would not allow for the scoring of user-provided signatures.
* Fixed the TB_hiv data to remove unnecessary factor level of Disease metadata. 
* Fixed the row numbers of existing sigAnnotData and common_sigAnnotData objects, and added code to update them after new signatures are added.

## Major Changes
* If any signatures in the object used with runTBsigProfiler() have <2 genes present in the given sample, the signatures will not be scored. This may affect existing scripts.

## Minor Changes
* Reorganized code for OriginalModel.R for clarity.
* Added 4 new signatures (Tabone_OD_11/TB12, Tabone_RES_25/EarlyRESP-TB25, Tabone_RES_27/TREAT-TB27, Long_RES_10)
* Updated the website interface
* Changed HGNChelper installation to be checked during profiling if update_genes = TRUE
* Reorganized code in mkAssay() for clarity. The output_name argument is now appended to all output assays, whereas previously it was only appended to the log of the input assay.

# TBSignatureProfiler 1.6.0

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