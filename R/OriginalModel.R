globalVariables(c("TBsignaturesSplit", "OriginalTrainingData", "TBsignatures"))

#' A function that implements the original methods for multiple TB signatures.
#'
#' This function computes prediction foot multiple TB signatures based on their training
#' models/methods. To avoid naming issues, the gene names for both training data and
#' input gene sets have been updated using the \code{\link[HGNChelper]{checkGeneSymbols}}.
#' TB signatures with available original models are: Anderson_42,
#' Anderson_OD_51, Kaforou_27, Kaforou_OD_44, Kaforou_OD_53, Sweeney_OD_3,
#' Maertzdorf_4, Verhagen_10, Jacobsen_3, Sambarey_HIV_10, Leong_24,
#' Berry_OD_86, Berry_393, Bloom_OD_144, Suliman_RISK_4, Zak_RISK_16,
#' Leong_RISK_29, and Zhao_NANO_6.
#' The predicted score for each signature has been stored in the column data
#' section of the input SummarizedExperiment study.
#'
#' @param input A SummarizedExperiment object with gene symbols as the assay row
#' names.
#' @param geneSignaturesName A character string/vector specifying the signature
#' of interest. If \code{any(geneSignaturesName == "") == TRUE} run all
#' available gene signatures' original models.
#' @param useAssay A character string or an integer specifying the assay in the
#' \code{input}. Default is the first assay in the assay list.
#' Used for the test SummarizedExperiment object. Default is \code{1}, indicating the
#' first assay in the \code{input}.
#' @param adj A small positive real number used in \code{\link[sva]{ComBat}} to solve
#' for genes with 0 counts(rare cases). Default is \code{1e-3}.
#' @param BPPARAM An instance inherited from \code{\link[BiocParallel]{bplapply}}.
#'
#' @return A SummarizedExperiment object with predicted scores for each sample
#' obtained from the signature's original model.
#' @export
#' @examples
#' data("TB_hiv", package = "TBSignatureProfiler")
#' re <- evaluateOriginalModel(input = TB_hiv,
#' geneSignaturesName = c("Anderson_42", "Maertzdorf_15"), useAssay = "counts")
#' re$Anderson_42_OriginalModel

evaluateOriginalModel <- function(input, geneSignaturesName, useAssay = 1, adj = 1e-3,
                                  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    ## Check input class, should be SummarizedExperiment
    if (!methods::is(input, "SummarizedExperiment")) {
        base::stop("The \"input\" should be a SummarizedExperiment object.")
    }
    signature_NoRetraining <- c("Anderson_42", "Anderson_OD_51", "Kaforou_27",
                                "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3")
    signature_Retraining <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                              "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                              "Leong_24", "Berry_OD_86", "Berry_393",
                              "Bloom_OD_144", "Suliman_RISK_4", "Zak_RISK_16",
                              "Leong_RISK_29", "Zhao_NANO_6")
    signature_all <- c(signature_NoRetraining, signature_Retraining)
    ## Check whether input geneSignaturesName is valid
    if (base::missing(geneSignaturesName)) {
        stop(base::sprintf("Argument \"geneSignaturesName\" is missing, with no default.\nAvailable signatures are: %s",
                           base::paste0(signature_all, collapse = ", ")))
    }
    if (base::any(geneSignaturesName == "")) {
        base::message("Evaluating all available signatures with their original model.")
        ## Remove "" from geneSignaturesName
        geneSignaturesName <- geneSignaturesName[-base::which(geneSignaturesName == "")] %>%
            base::union(signature_all)
    }
    ## Intersect input geneSignaturesName with available gene signatures
    sig_sub <- base::intersect(geneSignaturesName, signature_all)
    if (base::identical(sig_sub, character(0))) {
        stop(base::sprintf("Original model(s) for \"geneSignaturesName\": %s is/are not available.\nAvailable signatures are: %s",
                           base::paste0(geneSignaturesName, collapse = ", "),
                           base::paste0(signature_all, collapse = ", ")), call. = FALSE)
    }
    ## Check geneSignaturesName that are not available
    sig_not_found <- geneSignaturesName[-base::which(geneSignaturesName %in% sig_sub)]
    if (!base::identical(sig_not_found, character(0))) {
        base::message(base::sprintf("Original model(s) for \"geneSignaturesName\": %s are/is not available.",
                                    base::paste0(sig_not_found, collapse = ", ")))
    }
    sig_sub_NoRetraining <- base::intersect(sig_sub, signature_NoRetraining)
    if (!base::identical(sig_sub_NoRetraining, character(0))) {
        sample_score_NoRetraining <- .OriginalModel_NoRetraining(input, useAssay,
                                                                 sig_sub_NoRetraining,
                                                                 BPPARAM)
    } else {
        sample_score_NoRetraining <- NULL
    }
    sig_sub_Retraining <- base::intersect(sig_sub, signature_Retraining)
    if (!identical(sig_sub_Retraining, character(0))) {
        sample_score_Retraining <- .OriginalModel_Retraining(input, useAssay,
                                                             sig_sub_Retraining,
                                                             adj, BPPARAM)
    } else {
        sample_score_Retraining <- NULL
    }
    col_info <- SummarizedExperiment::colData(input)
    if (!base::is.null(sample_score_NoRetraining) && !base::is.null(sample_score_Retraining)) {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_NoRetraining,
                                                                 sample_score_Retraining)
    } else if (!base::is.null(sample_score_NoRetraining) && base::is.null(sample_score_Retraining)) {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_NoRetraining)
    } else {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_Retraining)
    }
    return(input)
}

#' TB gene signatures that do not require retraining.
#'
#' A function to obtain predicted score for TB gene signatures that do not need to be retrained.
#'
#' Anderson_42 and Anderson_OD_51 used difference of sums to calculate
#' prediction scores. Difference of sums is obtained by subtracting the sum of
#' the expression of genes within signatures that are down-regulated from the
#' sum of the expression of genes that are up-regulated within signatures.
#' Kaforou_27, Kaforou_OD_44, and Kaforou_OD_53 used difference of arithmetic
#' means to calculate prediction scores.
#' Sweeney_OD_3 used difference of arithmetic mean to calculate prediction score.
#'
#' @inheritParams evaluateOriginalModel
#' @return A SummarizedExperiment object with predicted scores for each sample
#' obtained from the signature's original model.

.OriginalModel_NoRetraining <- function(input, useAssay, geneSignaturesName, BPPARAM) {
    sample_score_list <- BiocParallel::bplapply(seq_len(length(geneSignaturesName)),
                                                function(i) {
                                                    TBsignaturesSplit <- TBsignaturesSplit
                                                    sig <- geneSignaturesName[i]
                                                    base::message(base::sprintf("Now Evaluating: %s", sig))
                                                    up_set <- TBsignaturesSplit[[sig]][[base::paste0(sig, "_up")]]
                                                    dn_set <- TBsignaturesSplit[[sig]][[base::paste0(sig, "_dn")]]
                                                    dat_sig_up <- base::data.frame(subsetGeneSet(input, up_set, useAssay))
                                                    dat_sig_dn <- base::data.frame(subsetGeneSet(input, dn_set, useAssay))
                                                    if (base::ncol(dat_sig_up) == 0 || base::ncol(dat_sig_dn) == 0) {
                                                        if (base::ncol(dat_sig_up) == 0) {
                                                            base::warning(base::sprintf("All up-regulated genes are not found within the input. Prediciton score for signature: %s is NA",
                                                                                        sig))
                                                        } else {
                                                            base::warning(base::sprintf("All down-regulated genes are not found within the input. Prediciton score for signature: %s is NA",
                                                                                        sig))
                                                        }
                                                        sample_score <- NA
                                                    } else { # Need number of column >=1
                                                        if (sig == "Anderson_42" || sig == "Anderson_OD_51") {
                                                            up_sample_score <- BiocGenerics::rowSums(dat_sig_up)
                                                            dn_sample_score <- BiocGenerics::rowSums(dat_sig_dn)
                                                        } else if (sig == "Kaforou_27" || sig == "Kaforou_OD_44"
                                                                   || sig == "Kaforou_OD_53") {
                                                            up_sample_score <- BiocGenerics::rowMeans(dat_sig_up)
                                                            dn_sample_score <- BiocGenerics::rowMeans(dat_sig_dn)
                                                        } else if (sig == "Sweeney_OD_3") {
                                                            ## Only one gene in down set for Sweeney_OD_3.
                                                            up_sample_score <- base::apply(dat_sig_up, 1, function(x) mean(x))
                                                            dn_sample_score <- base::apply(dat_sig_dn, 1, function(x) mean(x))
                                                        }
                                                        sample_score <- up_sample_score - dn_sample_score
                                                    }
                                                    sample_score
                                                }, BPPARAM = BPPARAM)
    sample_score_result <- base::do.call(base::cbind, sample_score_list)
    base::colnames(sample_score_result) <- base::paste0(geneSignaturesName, "_OriginalModel")
    return(sample_score_result)
}

#' TB gene signatures that require retraining.
#'
#' A function to obtain predicted score for TB gene signatures that need retraining of original models.
#'
#' Maertzdorf_4 and Maertzdorf_15 were trained using random forest to distinguish
#' patients with active TB from healthy controls.
#' Verhagen_10 was also trained using random forest to distinguish samples with active TB
#' from either latent infection or healthy controls.
#' The random forest model was build using \code{\link[randomForest]{randomForest}}.\cr
#' Jacobsen_3 were trained using linear discriminatory analysis (LDA)
#' to distinguish samples with active TB from latent infection status.
#' Sambarey_HIV_10 were also trained using LDA to distinguish samples with active TB
#' from either latent infection, healthy control, or other disease (HIV).
#' The LDA model was build using \code{\link[MASS]{lda}}.\cr
#' Berry_OD_86 and Berry_393 were trained using K-nearest neighbors (KNN) model to
#' differentiate samples with active TB from latent infection status.
#' The KNN model was build using \code{\link[class]{knn}}.\cr
#' Suliman_RISK_4 and Zak_RISK_16 were trained using support vector machines (SVM)
#' to distinguish TB progressor from non-progressors. The input gene expression features
#' for Suliman_RISK_4 used the paired ratio of GAS6/CD1C, SEPTIN4/BLK, SEPTIN4/CD1C, GAS6/BLK.
#' The SVM model was build using \code{\link[e1071]{svm}}.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @inheritParams evaluateOriginalModel
#' @return A SummarizedExperiment object with predicted scores for each sample
#' obtained from the signature's original model.

.OriginalModel_Retraining <- function(input, useAssay, geneSignaturesName, adj, BPPARAM) {
    OriginalTrainingData <- OriginalTrainingData
    TBsignatures <- TBsignatures
    gene_set_length <- base::length(geneSignaturesName)
    sample_score_list <- BiocParallel::bplapply(base::seq_len(gene_set_length), function(i) {
        ## Identify signature name
        sig <- geneSignaturesName[i]
        base::message(base::sprintf("Now Evaluating: %s", sig))
        if (sig == "Leong_RISK_29" || sig == "Zak_RISK_16") {
            theObject_train <- OriginalTrainingData[["GSE79362"]]
            if (sig == "Leong_RISK_29") {
                col_info <- SummarizedExperiment::colData(theObject_train)
                ## Only include baseline samples from GSE79362
                sample_baseline_GSE79362 <- col_info[, c("PatientID", "MeasurementTime")] %>%
                    base::as.data.frame() %>%
                    dplyr::mutate(sample_name = base::row.names(col_info)) %>%
                    dplyr::arrange(.data$MeasurementTime, .data$PatientID) %>%
                    dplyr::group_by(.data$PatientID) %>%
                    dplyr::mutate(first = dplyr::first(.data$sample_name))
                theObject_train <- theObject_train[, base::unique(sample_baseline_GSE79362$first)]
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input, SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "cv_glmnet_OriginalModel",
                                                                adj = adj)
            }
            if (sig == "Zak_RISK_16") {
                theObject_train <- theObject_train[, theObject_train$ACS_cohort %in% "Training"]
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "svm_OriginalModel",
                                                                adj = adj)
            }
            return(sample_score)
        } else if (sig == "Leong_24") {
            theObject_train <- OriginalTrainingData[["GSE101705"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "cv_glmnet_OriginalModel",
                                                            adj = adj)
            return(sample_score)
        } else if (sig == "Maertzdorf_4" || sig == "Maertzdorf_15") {
            theObject_train <- OriginalTrainingData[["GSE74092"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "randomForest_OriginalModel",
                                                            adj = adj)
            return(sample_score)
        }
        else if (sig == "LauxdaCosta_OD_3" || sig == "Bloom_OD_144") {
            theObject_train <- OriginalTrainingData[["GSE42834"]]
            if (sig == "LauxdaCosta_OD_3") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "randomForest_OriginalModel",
                                                                adj = adj)
            }
            if (sig == "Bloom_OD_144") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "svm_OriginalModel",
                                                                adj = adj)
            }
            return(sample_score)
        } else if (sig == "Jacobsen_3" || sig == "Berry_393" || sig == "Berry_OD_86") {
            theObject_train <- OriginalTrainingData[["GSE19491"]]
            if (sig == "Jacobsen_3") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "lda_OriginalModel",
                                                                adj = adj)
            }
            if (sig == "Berry_393" || sig == "Berry_OD_86") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "knn_OriginalModel",
                                                                adj = adj)
            }
            return(sample_score)
        } else if (sig == "Sambarey_HIV_10") {
            theObject_train <- OriginalTrainingData[["GSE37250"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "lda_OriginalModel",
                                                            adj = adj)
            return(sample_score)
        } else if (sig == "Verhagen_10") {
            theObject_train <- OriginalTrainingData[["GSE41055"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "randomForest_OriginalModel",
                                                            adj = adj)
            return(sample_score)
        } else if (sig == "Suliman_RISK_4") {
            theObject_train <- OriginalTrainingData[["GSE94438"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "SulimanOriginalModel",
                                                            adj = adj)
            return(sample_score)
        } else if (sig == "Zhao_NANO_6") {
            theObject_train <- OriginalTrainingData[["Zhao_NANO_6"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "cv_glmnet_OriginalModel",
                                                            adj = adj)
            return(sample_score)
        }
    }, BPPARAM = BPPARAM)
    sample_score_result <- base::do.call(base::cbind, sample_score_list)
    base::colnames(sample_score_result) <- base::paste0(geneSignaturesName, "_OriginalModel")
    return(sample_score_result)
}

#' Obtain training data, testing data, and train signature's original model.
#'
#' @inheritParams ref_combat_impute
#' @param obtainDiagnosis Boolean. Used to create training data if TRUE. Default is FALSE
#' @param annotationColName A character string specifying the column name of disease status.
#' Only used when creating training data. Default is NULL.
#' @param FUN A character string specifying the function name of the corresponding signature's original model.
#' @return The predicted score for each sample in the test study using corresponding gene signature's original model.
ObtainSampleScore_OriginalModel <- function(theObject_train, useAssay, gene_set,
                                            input, SigName, obtainDiagnosis,
                                            annotationColName, FUN, adj) {
    dat_test_sig <- ref_combat_impute(theObject_train, useAssay, gene_set, input,
                                      SigName, adj)
    if (base::is.null(dat_test_sig)) {
        sample_score <- base::rep(NA, base::ncol(input))
        return(sample_score)
    }
    ## Specify the useAssay = 1 argument for training data explicitly. DO NOT DELETE!
    dat_list <- subsetGeneSet(theObject_train, gene_set, useAssay = 1,
                              obtainDiagnosis, annotationColName)
    FUN <- base::match.fun(FUN)
    sample_score <- FUN(dat_list, dat_test_sig)
    return(sample_score)
}

#' Filter gene expression value matrix based on certain gene sets.
#'
#' A function used to subset gene expression value matrix based on certain gene sets.
#'
#' @param theObject A SummarizedExperiment object that has been prestored in the
#' RData file: OriginalTrainingData.
#' @param gene_set A character vector that includes gene symbols for gene signatures.
#' @param useAssay A character string or an integer specifying the assay in the \code{theObject}
#' that will be selected.
#' @param obtainDiagnosis Boolean. Usually used to create training data if TRUE. Default is FALSE
#' @param annotationColName A character string specifying the column name of disease status.
#' Only used when creating training data. Default is NULL.
#' @return A \code{matrix} with selected gene expression value if \code{obtainDiagnosis == FALSE}.
#' If \code{obtainDiagnosis == TRUE}, return a \code{list} contains the selected
#' gene expression value and diagnosis results for each sample.
#' @importFrom stats median na.pass
subsetGeneSet <- function(theObject, gene_set, useAssay,
                          obtainDiagnosis = FALSE, annotationColName = NULL) {
    dat_assay <- SummarizedExperiment::assays(theObject)[[useAssay]]
    ## Check assay names not found case
    if (base::is.null(dat_assay)) {
        stop(base::sprintf("Input object does not have assay with name: %s", useAssay))
    }
    dat_assay <- dat_assay %>%
        base::as.data.frame()
    ## Update gene names for both gene_set and row names of the input data
    ## Function update_genenames() were found from profile.R
    gene_set <- gene_set %>%
        update_genenames() %>%
        base::unique()
    update_rownames <- base::row.names(dat_assay) %>%
        update_genenames()
    ## Check for duplicated names after updating
    index <- base::which(update_rownames %in% gene_set)
    dat_assay_sig <- dat_assay[index, , drop = FALSE]
    if (base::max(base::table(update_rownames[index])) > 1) {
        ## Collapse duplicated gene names based on median value
        dat_assay_sig$SYMBOL <- update_rownames[index]
        exprs2 <- stats::aggregate(stats::as.formula(". ~ SYMBOL"),
                                   data = dat_assay_sig, FUN = median,
                                   na.action = na.pass)
        base::row.names(exprs2) <- exprs2$SYMBOL
        dat_assay_sig <- exprs2[, -base::which(base::colnames(exprs2) == "SYMBOL"),
                                drop = FALSE]
    } else {
        base::row.names(dat_assay_sig) <- update_rownames[index]
    }
    dat_assay_sig <- base::t(dat_assay_sig)
    dat_assay_sig <- dat_assay_sig[, base::sort(base::colnames(dat_assay_sig)),
                                   drop = FALSE] %>%
        base::as.matrix()
    if (obtainDiagnosis == FALSE) {
        return(dat_assay_sig)
    } else {
        if (base::is.null(annotationColName)) {
            stop("\"annotationColName\" is not given.
                 Need annotation column name to obtain sample diagnosis results.")
        }
        diagnosis_train <- SummarizedExperiment::colData(theObject)[, annotationColName]
        dat_list <- base::list(dat_train_sig = dat_assay_sig,
                               diagnosis_train = diagnosis_train)
        return(dat_list)
    }
}

#' A function for reference batch correction and imputation.
#'
#' A function used to perform reference batch correction and imputation in the
#' testing data for gene signatures that require retraining of the model.
#' We used the k-nearest neighbors to impute the expression values for missing gene(s).
#' The imputation operation is achieved using \code{\link[impute]{impute.knn}}.
#' Since the computational time for the imputation step can be excessive for large
#' number of missing genes. We made some constrains to prevent the overflow of imputation
#' operation. The evaluation will not run if more than \code{geneMax}*100\% percent
#' of the genes are not found for the corresponding gene signature in the input study.
#' By default \code{geneMax} = 0.8, so the evaluation will not run if more than 80\%
#' of the genes are missing when matching the input study to the reference data.
#'
#' @param theObject_train A SummarizedExperiment object that has been prestored in the
#' data file: OriginalTrainingData.
#' @param useAssay A character string or an integer specifying the assay in the \code{input}.
#' Used for the test SummarizedExperiment object. Default is 1, indicating the first assay
#' in the test SummarizedExperiment object.
#' @param gene_set A character vector that includes gene symbols for selected gene signature.
#' @param input A SummarizedExperiment object with gene symbols as the assay row names.
#' @param SigName Optional. A character string that indicates the name for \code{gene_set}.
#' \code{SigName} is used to provide information when gene signatures were missing
#' in the test data.
#' @param adj A small real number used in combat to solve for genes with 0 counts in rare cases.
#' Not required for most of cases.
#' @param geneMax A real number between 0 and 1. This is used to detect the
#' maximum percent missing genes allowed in the evaluated signatures.
#' See \code{\link[impute]{impute.knn}} for details. The default value is 0.8.
#' @return Gene set subset

ref_combat_impute <- function(theObject_train, useAssay, gene_set, input,
                              SigName, adj, geneMax = 0.8) {
    object_list <- base::list(GSE_train = theObject_train, GSE_test = input)
    ## Specify the useAssay = 1 argument for training data explicitly. DO NOT DELETE!
    dat_train_assay_sig <- base::t(subsetGeneSet(theObject_train, gene_set, 1))
    dat_test_assay_sig <- base::t(subsetGeneSet(input, gene_set, useAssay))
    index_match <- base::match(base::row.names(dat_train_assay_sig),
                               base::row.names(dat_test_assay_sig))
    index_match_boolean <- base::is.na(index_match)
    if (base::all(index_match_boolean)) {
        base::warning(base::sprintf("No genes are found in the input study within given gene signatures: %s. The model will not run and signature score will be NA",
                                    SigName))
        return(NULL)
    }
    if (base::sum(index_match_boolean) > base::floor(geneMax * base::length(index_match))) {
        base::warning(base::sprintf("More than %1.2f%% of the genes are missing in the input study \n for gene signatures: %s. ",
                                    100 * geneMax, SigName))
    }
    dat_exprs_match_sub <- base::list(GSE_train = base::data.frame(dat_train_assay_sig),
                                      GSE_test = base::data.frame(dat_test_assay_sig))
    dat_exprs_combine <- base::Reduce(function(x, y)
        base::merge(x, y, by = "id", all.x = TRUE),
        base::lapply(dat_exprs_match_sub, function(x) {
            x$id <- base::row.names(x)
            x
        }))
    row_names <- dat_exprs_combine$id
    dat_exprs_count <- dat_exprs_combine %>%
        dplyr::select(-.data$id) %>%
        base::as.data.frame()
    row.names(dat_exprs_count) <- row_names
    # Check for NA's in the dat_exprs_counts
    if (base::sum(base::is.na(dat_exprs_count)) == 0) {
        dat_exprs_count1 <- base::as.matrix(dat_exprs_count)
    } else {
        # Genes in the training data are not found in testing.
        # Use KNN imputation
        # DO some checks on the number of missing genes
        num_missing_gene <- base::sum(base::is.na(dat_exprs_count[, ncol(dat_exprs_count)]))
        # Calculate impute parameters for function: impute.knn
        colmax <- base::min(1, base::round(num_missing_gene / base::nrow(dat_exprs_count), 5) + 1e-5)
        rowmax <- base::min(1, base::round(ncol(dat_test_assay_sig) / base::ncol(dat_exprs_count), 5) + 1e-5)
        common_genes <- base::row.names(dat_test_assay_sig)[stats::na.omit(index_match)]
        missing_genes <- base::row.names(dat_train_assay_sig)[!base::row.names(dat_train_assay_sig) %in% common_genes]
        message(sprintf("Gene(s): %s is/are not found in the input study for gene signature: %s.\nImpute missing features using KNN",
                        paste0(missing_genes, collapse = ", "), SigName))
        output <- impute::impute.knn(as.matrix(dat_exprs_count), k = 10,
                                     rowmax = rowmax, colmax = colmax)
        dat_exprs_count1 <- output$data
    }
    ## Obtain annotation data from the samples
    col_data <- lapply(base::seq_len(base::length(object_list)), function(x) {
        col_data <- SummarizedExperiment::colData(object_list[[x]])
        col_data$GSE <- base::names(object_list[x])
        col_data
    })
    col_info <- plyr::rbind.fill(base::lapply(col_data, function(x) {
        x$SampleName <- base::row.names(x)
        base::as.data.frame(x)
    }))
    mod1 <- stats::model.matrix(~ 1, data = col_info)
    # Use reference combat with tryCatch
    # +1e-3 to solve for genes with 0 counts in rare cases
    ref_combat_tryCatch <- function(dat_exprs_count1) {
        tryCatch(
            expr = {
                sva::ComBat(dat = dat_exprs_count1,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
            },
            error = function(e) {
                sva::ComBat(dat = dat_exprs_count1 + adj,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
            }
        )
    }
    combat_edata <- ref_combat_tryCatch(dat_exprs_count1)
    runindata_test <- combat_edata[, (base::ncol(theObject_train) + 1):base::ncol(combat_edata)]
    base::colnames(runindata_test) <- base::colnames(input)
    ## Specify useAssay = 1 explicitly as we are creating new Sobject
    dat_test_sig <- SummarizedExperiment::SummarizedExperiment(assay = as.matrix(runindata_test)) %>%
        subsetGeneSet(gene_set, useAssay = 1)
    return(dat_test_sig)
}

#' Train original model for gene signatures Leong_24, Leong_RISK_29, Zhao_NANO_6 using lasso logistic regression.
#'
#' @param dat_list A \code{list} contains training data and disease status outcomes from
#' the discovery data of corresponding gene signatures.
#' @param dat_test_sig A \code{data frame} contains corresponding gene sets from the \code{input}.
#' @return The predicted score for each sample in the test study.

cv_glmnet_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB")) %>%
        base::as.integer()
    sig_model <- glmnet::cv.glmnet(x = dat_list$dat_train_sig,
                                   y = diagnosis_train)
    sample_score <- stats::predict(sig_model$glmnet.fit, s = sig_model$lambda.min,
                                   newx = dat_test_sig, na.action = na.omit)
    pred_score <- base::as.vector(sample_score)
    return(pred_score)
}

#' Train original model for gene signatures Maertzdorf_4, Maertzdorf_15, Verhagen_10,
#' and LauxdaCosta_OD_3.
#'
#' @inheritParams cv_glmnet_OriginalModel
#' @importFrom stats na.omit
#' @return The predicted score for each sample in the test study.

randomForest_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB"))
    base::set.seed(1)
    sig_model <- randomForest::randomForest(x = dat_list$dat_train_sig,
                                            y = as.factor(diagnosis_train),
                                            ntree = 5000, importance = TRUE)
    sample_score <- stats::predict(object = sig_model, newdata = dat_test_sig,
                                   type = "prob")
    pred_score <- base::unlist(sample_score[, 2], use.names = FALSE)
    return(pred_score)
}

#' Train original model for gene signatures Jacobsen_3 and Sambarey_HIV_10.
#'
#' @inheritParams cv_glmnet_OriginalModel
#' @return The predicted score for each sample in the test study.

lda_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB")) %>%
        base::as.integer()
    sig_model <- MASS::lda(diagnosis_train ~ .,
                           base::data.frame(dat_list$dat_train_sig, diagnosis_train))
    sample_score <- stats::predict(sig_model, base::data.frame(dat_test_sig))
    pred_score <- base::as.vector(sample_score$posterior[, 2])
    return(pred_score)
}

#' Train original model for gene signatures Berry_393 and Berry_OD_86.
#'
#' @inheritParams cv_glmnet_OriginalModel
#' @return The predicted score for each sample in the test study.

knn_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    dat_train_sig <- dat_list$dat_train_sig
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB")) %>%
        base::as.integer()
    sig_model <- class::knn(train = dat_train_sig, test = dat_test_sig,
                            cl = diagnosis_train, k = 10, prob = TRUE)
    sample_score_ori <- base::attributes(sig_model)$prob
    pred_result <- base::data.frame(sample_score_ori, sig_model)
    pred_result_nrow <- base::nrow(pred_result)
    pred_score <- rep(0, pred_result_nrow)
    for (i in base::seq_len(pred_result_nrow)) {
        pred_score[i] <- ifelse(pred_result$sig_model[i] == 2,
                                pred_result$sample_score_ori[i],
                                1 - pred_result$sample_score_ori[i])
    }
    return(pred_score)
}

#' Train original model for gene signatures Bloom_OD_144 and Zak_RISK_16.
#'
#' @inheritParams cv_glmnet_OriginalModel
#' @return The predicted score for each sample in the test study.

svm_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB")) %>%
        base::as.integer()
    sig_model <- suppressWarnings(e1071::svm(x = dat_list$dat_train_sig,
                                             y = diagnosis_train,
                                             type = "nu-regression", kernel = "linear",
                                             cost = 100, cachesize = 5000, tolerance = 0.01,
                                             shrinking = FALSE, cross = 3))
    pred_score <- stats::predict(sig_model, dat_test_sig)
    return(pred_score)
}

#' Train original model gene signature Suliman_RISK_4.
#'
#' @inheritParams cv_glmnet_OriginalModel
#' @return The predicted score for each sample in the test study.

SulimanOriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- base::factor(diagnosis_train_new, levels = c("Others", "PTB")) %>%
        base::as.integer()
    dat_train_sig <- dat_list$dat_train_sig
    dat_train_sig <- base::data.frame(GAS6_CD1C = dat_train_sig[, "GAS6"] / dat_train_sig[, "CD1C"],
                                      SEPTIN4_BLK = dat_train_sig[, "SEPTIN4"] / dat_train_sig[, "BLK"],
                                      SEPTIN4_CD1C = dat_train_sig[, "SEPTIN4"] / dat_train_sig[, "CD1C"],
                                      GAS6_BLK = dat_train_sig[, "GAS6"] / dat_train_sig[, "BLK"])
    dat_test_sig <- base::data.frame(GAS6_CD1C = dat_test_sig[, "GAS6"] / dat_test_sig[, "CD1C"],
                                     SEPTIN4_BLK = dat_test_sig[, "SEPTIN4"] / dat_test_sig[, "BLK"],
                                     SEPTIN4_CD1C = dat_test_sig[, "SEPTIN4"] / dat_test_sig[, "CD1C"],
                                     GAS6_BLK = dat_test_sig[, "GAS6"] / dat_test_sig[, "BLK"])
    sig_model <- suppressWarnings(e1071::svm(x = dat_train_sig,
                                             y = diagnosis_train,
                                             type = "nu-regression", kernel = "linear",
                                             cost = 100, cachesize = 5000, tolerance = 0.01,
                                             shrinking = FALSE, cross = 3))
    pred_score <- stats::predict(sig_model, dat_test_sig)
    return(pred_score)
}
