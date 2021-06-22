# TBSignatureProfiler (development version)

* Currently have 54 signatures available

## Bug Fixes
* Fixed incorrect names of signatures in `Tbcommon` and `common_sigAnnotData` objects

## Major Changes
* Added Estevez_133, Estevez_259, LauxdaCosta_OD_3, and Maertzdorf_15 signatures to the package
* Added Chen_HIV_4, Gliddon_HIV_3, Gliddon_2_OD_4, Kulkarni_HIV_2, and Heycken_FAIL_22 signatures
* Added a COVIDsignatures object to the package that can be used to profile COVID-19 gene transcript signatures, thanks to collaborator Dylan Sheerin (WEHI)
* Added functions to evaluate some signatures using their original models from Johnson lab member Xutao Wang

## Minor Changes
* Added `addTBsignature()` to more easily facilitate updating signatures in package
* Added pROC option to obtain confidence intervals on AUC values as part of \code{tableAUC()}
* Added citation for newly published paper

# TBSignatureProfiler 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.