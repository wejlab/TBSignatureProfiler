globalVariables(c("TBsignaturesSplit", "OriginalTrainingData"))

#' A function that implements TB gene signatures' original model.
#'
#' TB signatures with available original models are: Anderson_42,
#' Anderson_OD_51, Kaforou_27, Kaforou_OD_44, Kaforou_OD_53, Sweeney_OD_3,
#' Maertzdorf_4, Verhagen_10, Jacobsen_3, Sambarey_HIV_10, Leong_24, Berry_OD_86,
#' Berry_393, Bloom_OD_144, Suliman_RISK_4, Zak_RISK_16, and Leong_RISK_29.
#' The predicted score for each signature has been stored in the column data section
#' of the input SummarizedExperiment study.
#'
#' @param input A SummarizedExperiment object with gene symbols as the assay row names.
#' @param useAssay A character string or an integer specifying the assay in the \code{input}.
#' Used for the test SummarizedExperiment object. Default is 1, indicating the
#' first assay in the \code{input}.
#' @param geneSignaturesName A character string/vector specifying the signature of interest.
#' Default is NULL. If geneSignaturesName is NULL, run all available gene signatures' original model in the package.
#' @param BPPARAM An instance inherited from \code{bplappy}.
#' See \code{\link[BiocParallel]{bplapply}} for details.
#'
#' @return A SummarizedExperiment object with predicted scores for each sample
#' obtained from the signature's original model.
#' @export
#' @examples
#' library(TBSignatureProfiler)
#' data("TB_hiv",package = "TBSignatureProfiler")
#' re <- evaluateOriginalModel(input = TB_hiv, useAssay = "counts",
#' geneSignaturesName = c("Anderson_42","Maertzdorf_15"))
#' re$Anderson_42_OriginalModel
evaluateOriginalModel <- function(input, useAssay = 1, geneSignaturesName = NULL,
                                   BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    signature_NoRetraining <- c("Anderson_42", "Anderson_OD_51", "Kaforou_27",
                                "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3")
    signature_Retraining <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                              "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                              "Leong_24", "Berry_OD_86", "Berry_393", "Bloom_OD_144",
                              "Suliman_RISK_4", "Zak_RISK_16", "Leong_RISK_29")
    signature_all <- c(signature_NoRetraining, signature_Retraining)
    if (is.null(geneSignaturesName)) {
       geneSignaturesName <- signature_all
    }
    sig_sub <- base::intersect(geneSignaturesName, signature_all)
    if (identical(sig_sub, character(0))) {
       stop(sprintf("No available original model(s) information for the input signature(s). Available signatures are: %s",
                     paste0(signature_all, collapse = ", ")))
    }
    sig_not_found <- geneSignaturesName[-which(geneSignaturesName %in% sig_sub)]
    if (!identical(sig_not_found, character(0))) {
        message(sprintf("The original model for signature(s): %s are/is not available.",
                        paste0(sig_not_found, collapse = ", ")))
    }
    sig_sub_NoRetraining <- base::intersect(sig_sub, signature_NoRetraining)
    sig_sub_Retraining <- base::intersect(sig_sub, signature_Retraining)
    if (!identical(sig_sub_NoRetraining, character(0))) {
        sample_score_NoRetraining <- .OriginalModel_NoRetraining(input, useAssay,
                                                                 sig_sub_NoRetraining,
                                                                 BPPARAM)
    } else {
      sample_score_NoRetraining <- NULL
    }
    if (!identical(sig_sub_Retraining, character(0))) {
        sample_score_Retraining <- .OriginalModel_Retraining(input, useAssay,
                                                             sig_sub_Retraining,
                                                             BPPARAM)
    } else {
      sample_score_Retraining <- NULL
    }
    col_info <- SummarizedExperiment::colData(input)
    if (!is.null(sample_score_NoRetraining) && !is.null(sample_score_Retraining)) {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_NoRetraining,
                                                                 sample_score_Retraining)
    } else if (!is.null(sample_score_NoRetraining) && is.null(sample_score_Retraining)) {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_NoRetraining)
    } else {
        SummarizedExperiment::colData(input) <- S4Vectors::cbind(col_info,
                                                                 sample_score_Retraining)
    }
    return(input)
}

#' TB gene signatures do not require retraining
#' A function to obtain predicted score for TB gene signatures that do not need to retrain.
#'
#' Anderson_42 and Anderson_OD_51 used difference of sums to calculate prediction scores.
#' Difference of sums is obtained by subtracting the sum of the expression of
#' genes within signatures that are down-regulated from the sum of the expression
#' of genes that are up-regulated within signatures.\cr
#' Kaforou_27, Kaforou_OD_44, and Kaforou_OD_53 used difference of arithmetic means
#' to calculate prediction scores.\cr
#' Sweeney_OD_3 used difference of geometric mean to calculate prediction score.
#'
#' @inheritParams evaluateOriginalModel
.OriginalModel_NoRetraining <- function(input, useAssay, geneSignaturesName, BPPARAM) {
    sample_score_list <- BiocParallel::bplapply(seq_len(length(geneSignaturesName)), function(i) {
        TBsignaturesSplit <- TBSignatureProfiler::TBsignaturesSplit
        sig <- geneSignaturesName[i]
        message(sprintf("Now Evaluating: %s", sig))
        up_set <- TBsignaturesSplit[[sig]][[paste0(sig, "_up")]]
        dn_set <- TBsignaturesSplit[[sig]][[paste0(sig, "_dn")]]
        dat_sig_up <- data.frame(subsetGeneSet(input, up_set, useAssay))
        dat_sig_dn <- data.frame(subsetGeneSet(input, dn_set, useAssay))
        if (ncol(dat_sig_up) == 0) {
            message("No signature genes found in ", paste(up_set, collapse = " "),
                    " for signature ", sig, ", prediction score is NA.")
            up_sample_score <- NA
            } else if (ncol(dat_sig_dn) == 0) {
            message("No signature genes found in ", paste(dn_set, collapse = " "),
                    " for signature ", sig, ", prediction score is NA.")
            dn_sample_score <- NA
            }
        else { # Need number of column >=1
            if (sig == "Anderson_42" || sig == "Anderson_OD_51") {
                up_sample_score <- BiocGenerics::rowSums(dat_sig_up)
                dn_sample_score <- BiocGenerics::rowSums(dat_sig_dn)
            } else if (sig == "Kaforou_27" || sig == "Kaforou_OD_44"
                    || sig == "Kaforou_OD_53") {
                up_sample_score <- BiocGenerics::rowMeans(dat_sig_up)
                dn_sample_score <- BiocGenerics::rowMeans(dat_sig_dn)
            } else if (sig == "Sweeney_OD_3") { # Only one gene in down set.
                up_sample_score <- base::apply(dat_sig_up, 1, function(x) exp(mean(log(x))))
                dn_sample_score <- base::apply(dat_sig_dn, 1, function(x) exp(mean(log(x))))
            }
        }
        sample_score <- up_sample_score - dn_sample_score
        sample_score
    }, BPPARAM = BPPARAM)
    sample_score_result <- do.call(cbind, sample_score_list)
    colnames(sample_score_result) <- paste0(geneSignaturesName, "_OriginalModel")
    return(sample_score_result)
}

#' TB gene signatures require retraining
#'
#' A function to obtain predicted score for TB gene signatures that required retraining of their original model.
#'
#' Maertzdorf_4 and Maertzdorf_15 were trained using random forest to distinguish
#' samples with active TB from healthy controls.
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
#' for Suliman_RISK_4 used the paired ratio of GAS6/CD1C, SEPT4/BLK, SEPT4/CD1C, GAS6/BLK.
#' The SVM model was build using \code{\link[e1071]{svm}}.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @inheritParams evaluateOriginalModel
.OriginalModel_Retraining <- function(input, useAssay, geneSignaturesName, BPPARAM) {
    OriginalTrainingData <- TBSignatureProfiler::OriginalTrainingData
    TBsignatures <- TBSignatureProfiler::TBsignatures
    sample_score_list <- BiocParallel::bplapply(seq_len(length(geneSignaturesName)), function(i) {
        sig <- geneSignaturesName[i]  # Identify signature name
        message(sprintf("Now Evaluating: %s", sig))
        if (sig == "Leong_RISK_29" || sig == "Zak_RISK_16") {
            theObject_train <- OriginalTrainingData[["GSE79362"]] # Only include baseline samples from GSE79362
            if (sig == "Leong_RISK_29") {
                col_info <- SummarizedExperiment::colData(theObject_train)
                index_NA <- which(is.na(col_info$PatientID))
                if (requireNamespace("dplyr", quietly = TRUE)) {
                  sample_baseline_GSE79362 <- col_info[-index_NA, c("PatientID", "MeasurementTime")] %>%
                    data.frame() %>%
                    dplyr::mutate(sample_name = row.names(col_info)[-index_NA]) %>%
                    dplyr::arrange(.data$MeasurementTime, .data$PatientID) %>%
                    dplyr::group_by(.data$PatientID) %>%
                    dplyr::mutate(first = dplyr::first(.data$sample_name))
                } else {
                  message("Installing packages: \"dplyr\" for the function to work")
                  utils::install.packages("dplyr")
                  sample_baseline_GSE79362 <- col_info[-index_NA, c("PatientID", "MeasurementTime")] %>%
                      data.frame() %>%
                      dplyr::mutate(sample_name = row.names(col_info)[-index_NA]) %>%
                      dplyr::arrange(.data$MeasurementTime, .data$PatientID) %>%
                      dplyr::group_by(.data$PatientID) %>%
                      dplyr::mutate(first = dplyr::first(.data$sample_name))
                  }
                theObject_train <- theObject_train[, unique(sample_baseline_GSE79362$first)]
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input, SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "LeongOriginalModel")
            }
            if (sig == "Zak_RISK_16") {
                theObject_train <- theObject_train[, theObject_train$ACS_cohort %in% "Training"]
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input, SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "Bloom_Zak_OriginalModel")
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
                                                            FUN = "LeongOriginalModel")
            return(sample_score)
        } else if (sig == "Maertzdorf_4" || sig == "Maertzdorf_15") {
            theObject_train <- OriginalTrainingData[["GSE74092"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "Maertzdorf_Verhagen_daCosta_OriginalModel")
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
                                                                FUN = "Maertzdorf_Verhagen_daCosta_OriginalModel")
            }
            if (sig == "Bloom_OD_144") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "Bloom_Zak_OriginalModel")
            }
            return(sample_score)
        }
        else if (sig == "Jacobsen_3" || sig == "Berry_393" || sig == "Berry_OD_86") {
            theObject_train <- OriginalTrainingData[["GSE19491"]]
            if (sig == "Jacobsen_3") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "Jacobsen_Sambarey_OriginalModel")
                }
            if (sig == "Berry_393" || sig == "Berry_OD_86") {
                sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                                gene_set = TBsignatures[[sig]],
                                                                input = input,
                                                                SigName = sig,
                                                                obtainDiagnosis = TRUE,
                                                                annotationColName = "TBStatus",
                                                                FUN = "BerryOriginalModel")
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
                                                            FUN = "Jacobsen_Sambarey_OriginalModel")
            return(sample_score)
        } else if (sig == "Verhagen_10") {
            theObject_train <- OriginalTrainingData[["GSE41055"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "Maertzdorf_Verhagen_daCosta_OriginalModel")
            return(sample_score)
        } else if (sig == "Suliman_RISK_4") {
            theObject_train <- OriginalTrainingData[["GSE94438"]]
            sample_score <- ObtainSampleScore_OriginalModel(theObject_train, useAssay,
                                                            gene_set = TBsignatures[[sig]],
                                                            input = input,
                                                            SigName = sig,
                                                            obtainDiagnosis = TRUE,
                                                            annotationColName = "TBStatus",
                                                            FUN = "SulimanOriginalModel")
            return(sample_score)
        }
    }, BPPARAM = BPPARAM)
    sample_score_result <- do.call(cbind, sample_score_list)
    colnames(sample_score_result) <- paste0(geneSignaturesName, "_OriginalModel")
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
                                            annotationColName, FUN) {
    dat_test_sig <- ref_combat_impute(theObject_train, useAssay, gene_set, input, SigName)
    if (is.null(dat_test_sig)) {
        sample_score <- rep(NA, ncol(input))
        return(sample_score)
    }
    dat_list <- subsetGeneSet(theObject_train, gene_set, useAssay = 1,
                              obtainDiagnosis, annotationColName)
    FUN <- base::match.fun(FUN)
    sample_score <- FUN(dat_list, dat_test_sig)
    return(sample_score)
}

#' Filter gene expression value matrix based on certatin gene sets.
#'
#' A function used to subset gene expression value matrix based on certatin gene sets.
#'
#' @param theObject A SummarizedExperiment object that has been prestored in the
#' RData file: OriginalTrainingData.
#' @param gene_set A character vector that includes gene symbols for gene signatures.
#' @param useAssay A character string or an integer specifying the assay in the \code{theObject}
#' that will be selected. Default is 1, indicating the first assay of the SummarizedExperiment object.
#' @param obtainDiagnosis Boolean. Usually used to create training data if TRUE. Default is FALSE
#' @param annotationColName A character string specifying the column name of disease status.
#' Only used when creating training data. Default is NULL.
#' @return A matrix with selected gene expression value if \code{obtainDiagnosis} is FALSE.
#' else, return a list contains the selected gene expression value and diagnosis results for each observation.
subsetGeneSet <- function(theObject, gene_set, useAssay = 1,
                          obtainDiagnosis = FALSE, annotationColName = NULL) {
    dat_assay <- base::data.frame(SummarizedExperiment::assay(theObject, useAssay))
    dat_assay_sig <- base::t(dat_assay[which(base::row.names(dat_assay) %in% gene_set), ])
    dat_assay_sig <- dat_assay_sig[, base::sort(colnames(dat_assay_sig))]
    if (obtainDiagnosis == FALSE) {
      return(dat_assay_sig)
      } else {
        if (is.null(annotationColName)) {
            stop("annotationColName is not given.
                 Need annotation column name to obtain sample diagnosis results.")
        }
        diagnosis_train <- SummarizedExperiment::colData(theObject)[, annotationColName]
        dat_list <- list(dat_train_sig = dat_assay_sig,
                         diagnosis_train = diagnosis_train)
        return(dat_list)
    }
}

#' A function for reference batch correction and imputation
#'
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
#' By default \code{geneMax}=0.8, so the evaluation will not run if more than 80\%
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
#' in the test data. The default is an empty string.
#' @param geneMax A real number between 0 and 1. This is used to detect the
#' maximum percent missing genes allowed in the evaluated signatures.
#' See \code{\link[impute]{impute.knn}} for details. The default value is 0.8.
ref_combat_impute <- function(theObject_train, useAssay = 1, gene_set, input,
                              SigName = "", geneMax = 0.8) {
    object_list <- list(GSE_train = theObject_train, GSE_test = input)
    dat_train_assay_sig <- base::t(subsetGeneSet(theObject_train, gene_set, useAssay = 1))
    dat_test_assay_sig <- base::t(subsetGeneSet(input, gene_set, useAssay))
    index_match <- base::match(row.names(dat_train_assay_sig),
                               row.names(dat_test_assay_sig))
    if (identical(stats::na.omit(index_match), logical(0))) {
        message(sprintf("No genes are found in the input study within given gene signatures: %s. The model will not run and signature score will be NA",
                        SigName))
        return(NULL)
    }
    if (sum(is.na(index_match)) > floor(geneMax * length(index_match))) {
        warning(sprintf("More than %1.2f%% of the genes are missing in the input study \n for gene signatures: %s. ",
                        100 * geneMax, SigName))
    }
    dat_exprs_match_sub <- list(GSE_train = data.frame(dat_train_assay_sig),
                                GSE_test = data.frame(dat_test_assay_sig))
    dat_exprs_combine <- Reduce(
        function(x, y) merge(x, y, all.x = TRUE, by = "id"),
        lapply(dat_exprs_match_sub, function(x) {
          x$id <- rownames(x)
          x })
        )
    row_names <- dat_exprs_combine$id
    dat_exprs_count <- data.frame(dat_exprs_combine[, -which(colnames(dat_exprs_combine) == "id")])
    row.names(dat_exprs_count) <- row_names
    
    # Check for NA's in the dat_exprs_counts
    if (sum(is.na(dat_exprs_count)) == 0) {
        dat_exprs_count1 <- as.matrix(dat_exprs_count)
    } else {
        # Genes in the training data are not found in testing.
        # Use KNN imputation
        # DO some checks on the number of missing genes
        num_missing_gene <- sum(is.na(dat_exprs_count[, ncol(dat_exprs_count)]))
        # Calculate impute parameters for function: impute.knn
        colmax <- min(1, round(num_missing_gene / nrow(dat_exprs_count), 5) + 1e-5)
        rowmax <- min(1, round(ncol(dat_test_assay_sig) / ncol(dat_exprs_count), 5) + 1e-5)
        common_genes <- row.names(dat_test_assay_sig)[na.omit(index_match)]
        missing_genes <- row.names(dat_train_assay_sig)[!row.names(dat_train_assay_sig)
                                                        %in% common_genes]
        message(sprintf("Gene(s): %s is(are) not found in the input study for gene signature: %s. \n Impute missing features using KNN",
                        paste0(missing_genes, collapse = ", "), SigName))
        
        if (requireNamespace("impute", quietly = TRUE)) {
          output <- impute::impute.knn(as.matrix(dat_exprs_count), k = 10,
                                       rowmax = rowmax, colmax = colmax)
        } else {
          message("Installing required packages: \"impute\" for the function to work")
          utils::installed.packages("impute")
          output <- impute::impute.knn(as.matrix(dat_exprs_count), k = 10,
                                       rowmax = rowmax, colmax = colmax)
        }
        dat_exprs_count1 <- output$data
    }
    col_data <- lapply(seq_len(length(object_list)), function(x) { # Obtain annotation data from the samples
        col_data <- SummarizedExperiment::colData(object_list[[x]])
        col_data$GSE <- names(object_list[x])
        col_data
    })
    if (requireNamespace("plyr", quietly = TRUE)) {
      col_info <- plyr::rbind.fill(lapply(col_data, function(x) {
        x$SampleName <- row.names(x)
        as.data.frame(x)
      }))
    } else {
      message("Installing required package: \"plyr\" for the function to work")
      utils::installed.packages("plyr")
      col_info <- plyr::rbind.fill(lapply(col_data, function(x) {
        x$SampleName <- row.names(x)
        as.data.frame(x)
      }))
    }
    mod1 <- stats::model.matrix(~ 1, data = col_info)
    # Use reference combat with tryCatch
    # +1e-3 to solve for genes with 0 counts in rare cases
    ref_combat_tryCatch <- function(dat_exprs_count1) {
        tryCatch(
            expr = {
              if (requireNamespace("sva", quietly = TRUE)) {
                sva::ComBat(dat = dat_exprs_count1,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
              } else {
                message("Installing required packages: \"sva\" for the function to work")
                utils::installed.packages("sva")
                sva::ComBat(dat = dat_exprs_count1,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
              }
            },
            error = function(e) {
              if (requireNamespace("sva", quietly = TRUE)) {
                sva::ComBat(dat = dat_exprs_count1 + 1e-3,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
              } else {
                message("Installing required packages: \"sva\" for the function to work")
                utils::installed.packages("sva")
                sva::ComBat(dat = dat_exprs_count1 + 1e-3,
                            batch = col_info$GSE, mod = mod1,
                            ref.batch = "GSE_train")
              }
            }
        )
    }
    combat_edata <- ref_combat_tryCatch(dat_exprs_count1)
    runindata_test <- combat_edata[, (ncol(theObject_train) + 1):ncol(combat_edata)]
    colnames(runindata_test) <- colnames(input)
    theObject_test <- SummarizedExperiment::SummarizedExperiment(assay = as.matrix(runindata_test))
    dat_test_sig <- subsetGeneSet(theObject_test, gene_set, useAssay = 1)
    return(dat_test_sig)
}

#' Train original model for gene signatures: Leong_24 and Leong_RISK_29
#'
#' @param dat_list A list contains training data and disease status outcomes from
#' the discovery data of corresponding gene signatures.
#' @param dat_test_sig A data frame contains correpsonding biomarker gene features.
#' Used to evaluate the performance of biomarker using its original models.
#' @return The predicted score for each sample in the test study.
LeongOriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- as.integer(factor(diagnosis_train_new,
                                         levels = c("Others", "PTB")))
    if (requireNamespace("glmnet", quietly = TRUE)) {
      sig_model <- glmnet::cv.glmnet(x = dat_list$dat_train_sig,
                                     y = diagnosis_train)
    } else {
      message("Installing required packages: \"glmnet\" for the function to work")
      utils::installed.packages("glmnet")
      sig_model <- glmnet::cv.glmnet(x = dat_list$dat_train_sig,
                                     y = diagnosis_train)
    }
    sample_score <- stats::predict(sig_model$glmnet.fit, s = sig_model$lambda.min,
                                   newx = dat_test_sig, na.action = na.omit)
    pred_score <- as.vector(sample_score)
    return(pred_score)
}

#' Train original model for gene signatures: Maertzdorf_4, Maertzdorf_15, Verhagen_10,
#' and LauxdaCosta_OD_3
#'
#' @inheritParams LeongOriginalModel
#' @importFrom stats na.omit
#' @return The predicted score for each sample in the test study.
Maertzdorf_Verhagen_daCosta_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- factor(diagnosis_train_new, levels = c("Others", "PTB"))
    base::set.seed(1)
    if (requireNamespace("randomForest", quietly = TRUE)) {
      sig_model <- randomForest::randomForest(x = dat_list$dat_train_sig,
                                              y = as.factor(diagnosis_train),
                                              ntree = 5000, importance = TRUE)
    } else {
      message("Installing required packages: \"randomForest\" for the function to work")
      utils::installed.packages("randomForest")
      sig_model <- randomForest::randomForest(x = dat_list$dat_train_sig,
                                              y = as.factor(diagnosis_train),
                                              ntree = 5000, importance = TRUE)
    }
    sample_score <- stats::predict(object = sig_model, newdata = dat_test_sig,
                                   type = "prob")
    pred_score <- unlist(sample_score[, 2], use.names = FALSE)
    return(pred_score)
}

#' Train original model for gene signatrues: Jacobsen_3 and Sambarey_HIV_10
#'
#' @inheritParams LeongOriginalModel
#' @return The predicted score for each sample in the test study.
Jacobsen_Sambarey_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- as.integer(factor(diagnosis_train_new, levels = c("Others", "PTB")))
    if (requireNamespace("MASS", quietly = TRUE)) {
      sig_model <- MASS::lda(diagnosis_train ~ ., data.frame(dat_list$dat_train_sig,
                                                             diagnosis_train))
    } else {
      message("Installing required packages: \"MASS\" for the function to work")
      utils::installed.packages("MASS")
      sig_model <- MASS::lda(diagnosis_train ~ ., data.frame(dat_list$dat_train_sig,
                                                             diagnosis_train))
    }
    sample_score <- stats::predict(sig_model, data.frame(dat_test_sig))
    pred_score <- as.vector(sample_score$posterior[, 2])
    return(pred_score)
}

#' Train original model for gene signatures: Berry_393 and Berry_OD_86
#'
#' @inheritParams LeongOriginalModel
#' @return The predicted score for each sample in the test study.
BerryOriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    dat_train_sig <- dat_list$dat_train_sig
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- as.integer(factor(diagnosis_train_new, levels = c("Others", "PTB")))
    if (requireNamespace("class", quietly = TRUE)) {
      sig_model <- class::knn(train = dat_train_sig, test = dat_test_sig,
                              cl = diagnosis_train, k = 10, prob = TRUE)
    } else {
      message("Installing required packages: \"class\" for the function to work")
      utils::installed.packages("class")
      sig_model <- class::knn(train = dat_train_sig, test = dat_test_sig,
                              cl = diagnosis_train, k = 10, prob = TRUE)
    }
    sample_score_ori <- attributes(sig_model)$prob
    pred_result <- data.frame(sample_score_ori, sig_model)
    pred_score <- c()
    for (i in seq_len(nrow(pred_result))) {
        if (pred_result$sig_model[i] == 2) {
            pred_score[i] <- pred_result$sample_score_ori[i]
        } else {
            pred_score[i] <- 1 - pred_result$sample_score_ori[i]
        }
    }
    return(pred_score)
}

#' Train original model for gene signatures Bloom_OD_144 and Zak_RISK_16
#'
#' @inheritParams LeongOriginalModel
#' @return The predicted score for each sample in the test study.
Bloom_Zak_OriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- as.integer(factor(diagnosis_train_new, levels = c("Others", "PTB")))
    if (requireNamespace("e1071", quietly = TRUE)) {
      sig_model <- suppressWarnings(e1071::svm(x = dat_list$dat_train_sig,
                                               y = diagnosis_train,
                                               type = "nu-regression", kernel = "linear",
                                               cost = 100, cachesize = 5000, tolerance = 0.01,
                                               shrinking = FALSE, cross = 3))
    } else {
      message("Installing required packages: \"e1071\" for the function to work")
      utils::installed.packages("e1071")
      sig_model <- suppressWarnings(e1071::svm(x = dat_list$dat_train_sig,
                                               y = diagnosis_train,
                                               type = "nu-regression", kernel = "linear",
                                               cost = 100, cachesize = 5000, tolerance = 0.01,
                                               shrinking = FALSE, cross = 3))
    }
    pred_score <- stats::predict(sig_model, dat_test_sig)
    return(pred_score)
}

#' Train original model gene signature Suliman_RISK_4
#'
#' @inheritParams LeongOriginalModel
#' @return The predicted score for each sample in the test study.
SulimanOriginalModel <- function(dat_list, dat_test_sig) {
    diagnosis_train <- dat_list$diagnosis_train
    diagnosis_train_new <- base::ifelse(diagnosis_train == "PTB", "PTB", "Others")
    diagnosis_train <- as.integer(factor(diagnosis_train_new, levels = c("Others", "PTB")))
    dat_train_sig <- dat_list$dat_train_sig
    dat_train_sig <- data.frame(GAS6_CD1C = dat_train_sig[, "GAS6"] / dat_train_sig[, "CD1C"],
                                SEPT4_BLK = dat_train_sig[, "SEPT4"] / dat_train_sig[, "BLK"],
                                SEPT4_CD1C = dat_train_sig[, "SEPT4"] / dat_train_sig[, "CD1C"],
                                GAS6_BLK = dat_train_sig[, "GAS6"] / dat_train_sig[, "BLK"])
    dat_test_sig <- data.frame(GAS6_CD1C = dat_test_sig[, "GAS6"] / dat_test_sig[, "CD1C"],
                               SEPT4_BLK = dat_test_sig[, "SEPT4"] / dat_test_sig[, "BLK"],
                               SEPT4_CD1C = dat_test_sig[, "SEPT4"] / dat_test_sig[, "CD1C"],
                               GAS6_BLK = dat_test_sig[, "GAS6"] / dat_test_sig[, "BLK"])
    if (requireNamespace("e1071", quietly = TRUE)) {
      sig_model <- suppressWarnings(e1071::svm(x = dat_train_sig,
                                               y = diagnosis_train,
                                               type = "nu-regression", kernel = "linear",
                                               cost = 100, cachesize = 5000, tolerance = 0.01,
                                               shrinking = FALSE, cross = 3))
    } else {
      message("Installing required packages: \"e1071\" for the function to work")
      utils::installed.packages("e1071")
      sig_model <- suppressWarnings(e1071::svm(x = dat_train_sig,
                                               y = diagnosis_train,
                                               type = "nu-regression", kernel = "linear",
                                               cost = 100, cachesize = 5000, tolerance = 0.01,
                                               shrinking = FALSE, cross = 3))
    }
    pred_score <- stats::predict(sig_model, dat_test_sig)
    return(pred_score)
}
