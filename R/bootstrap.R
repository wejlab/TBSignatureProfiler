#' Bootstrap the AUC and Conduct T-Tests for a Collection of Signatures
#'
#' Run bootstrapping of the AUC and derive the p-value for a 2-sample t-test
#' for all signatures tested on a given dataset.
#' 
#' @param SE_scored a \code{SummarizedExperiment} object with genes as the 
#' row features and signature scores in the \code{colData}. There should also 
#' be a column of annotation data. Required.
#' @param annotationColName a character string giving the column name in 
#' \code{colData} that contains the annotation data. Required.
#' @param SignatureColNames a vector of column names in the 
#' \code{colData} that contain the signature score data. Required.
#' @param num.boot integer. The number of times to bootstrap the data. The 
#' default is \code{100}.
#' @param pb.show logical for whether to show a progress bar while running code.
#' The default is \code{TRUE}.
#' 
#' @return A list of length 3 returning a vector of p-values for a 2-sample 
#' t-test, bootstrapped AUC values, and an AUC value for using all scored values 
#' for all signatures specified in \code{SignatureColNames}.
#' 
#' @export
#' 
#' @examples 
#'  # Run signature profiling
#'  choose_sigs <- TBsignatures[-c(18, 22)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",  
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs)
#'  # Bootstrapping
#'  booted <- bootstrapAUC(SE_scored = prof_indian, annotationColName = "label",
#'                         signatureColNames = names(choose_sigs))
#'  str(booted)
#'  
bootstrapAUC <- function(SE_scored, annotationColName, signatureColNames, 
                         num.boot = 100, pb.show = TRUE){
  pvals <- aucs <- aucs_boot <- NULL
  
  annotationData <- colData(SE_scored)[annotationColName][, 1]
  if (!is.factor(annotationData)) annotationData <- as.factor(annotationData)
  
  # Create progress bar
  total <- length(signatureColNames)
  counter <- 0
  if(pb.show)  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in signatureColNames) {
    score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
    # Conduct a 2-sample t-test on the scores and their
    # corresponding Tuberculosis group status
    pvals <- c(pvals, t.test(score ~ annotationData)$p.value)
    
    # Obtain AUC based on entire dataset
    pred <- ROCit::rocit(score, annotationData)
    auc <- pred$AUC
    aucs <- c(aucs, max(auc, 1 - auc))
    
    # Proceed with bootstrapping
    tmp_aucs <- NULL
    for (j in 1:num.boot) {
      index <- sample(1:length(score), replace = T)
      tmp_score <- score[index]
      tmp_annotationData <- annotationData[index]
      pred <- ROCit::rocit(tmp_score, tmp_annotationData)
      tmp_auc <- max(pred$AUC, 1 - pred$AUC)
      tmp_aucs <- c(tmp_aucs, tmp_auc)
    }
    aucs_boot <- cbind(aucs_boot, tmp_aucs)
    
    # Update the progress bar
    counter <- counter + 1
    if (pb.show) setTxtProgressBar(pb, counter)
  }
  
  if (pb.show) close(pb)
  return(list("P-values" = pvals, "Boot AUC Values" = aucs_boot,
         "Non-Boot AUC Values" = aucs))
}

#' Create a Table of Results for T-tests and Rootstrapped AUC for a Collection of Scored Signatures.
#' 
#' This function collects the results of bootstrapping and t-tests for a scored 
#' gene expression dataset and presents them using a JavaScript table with an 
#' R interface.
#'
#' @inheritParams bootstrapAUC
#' 
#' @export
#' 
#' @return A JavaScript table with an R interface using the \code{DT} package. 
#' 
#' @examples 
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[-c(18, 22)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",  
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs)
#'  # Create table
#'  tableAUC(SE_scored = prof_indian, annotationColName = "label", 
#'           signatureColNames = names(choose_sigs))
#'  
tableAUC <- function(SE_scored, annotationColName, signatureColNames, 
                     num.boot = 100, pb.show = TRUE){
  # Run the bootstrapping function
  BS.Results <- bootstrapAUC(SE_scored, annotationColName, signatureColNames, 
                             num.boot, pb.show)
  pvals <- BS.Results[["P-values"]]
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  
  # Create interactive table
  return(DT::datatable(cbind("Signatures" = signatureColNames, 
                             "P-values" = round(pvals, 4), 
                             "-10*Log(p-values)" = round(-10 * log(pvals), 4), 
                             "LowerAUC" = round(apply(aucs_boot, 2, quantile, 
                                                      probs = .05), 4), 
                             "AUCs" = round(aucs, 4), 
                             "UpperAUC"= round(apply(aucs_boot, 2, quantile, 
                                                     probs = .95), 4)
                             )[order(aucs, decreasing = T),], 
                       options = list(scrollX = T, pageLength = 10), 
                       rownames = F))
}

#' Create a Comparison Plot of Boxplots for Bootstrapped AUC Values.
#' 
#' Present the results of AUC bootstrapping for a collection of scored 
#' signatures via boxplots.
#'
#' @inheritParams bootstrapAUC
#' @param num.boot an integer indicating the number of times to bootstrap the 
#' data.
#' @param cex.axis a number between 0 and 1 specifying the magnification to be 
#' used for axis annotation relative to the current setting of \code{cex}. 
#' @param cex a number between 0 and 1 giving the amount by which plotting text 
#' and symbols should be magnified relative to the default.
#' @param name a character string giving the overall title for the plot.
#' @param pb.show logical for whether to show a progress bar while running code.
#' 
#' @export
#' 
#' @return A plot with side-by-side boxplots of bootstrapped AUC values for 
#' each specified signature.
#' 
#' @examples 
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[-c(18, 22)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",  
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs)
#'  # Create boxplots
#'  compareBoxplots(prof_indian, "label", names(choose_sigs),
#'                  cex.axis = 0.5)
#'           
compareBoxplots <- function(SE_scored, annotationColName, signatureColNames, 
                            num.boot = 100, cex.axis = 0.7, cex = 0.25,
                            name = "Boxplot Comparison of Signatures",
                            pb.show = TRUE){
  BS.Results <- bootstrapAUC(SE_scored, annotationColName, signatureColNames, 
                             num.boot, pb.show)
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  colnames(aucs_boot) <- signatureColNames
  boxplot(aucs_boot[, order(aucs)], las = 2, cex.axis = cex.axis, cex = cex,
          main = name, ylab = "Bootstrapped AUC Values")
  abline(h = 0.5, col = "red", lty = 2)
}

#' Create an Array of ROC Plots to Compare Signatures.
#' 
#' @inheritParams signatureBoxplot
#' @param choose_colors a vector of length 2 definining the colors to be used
#' in the ROC plots. The default is c("cornflowerblue", "grey24").
#' 
#' @return An array of ROC plots.
#' 
#' @export
#' 
#' @examples 
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[-c(18, 22)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",  
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs)
#'                                  
#' # Create ROC plots
#' signatureROCplot(prof_indian, signatureColNames = names(choose_sigs),
#'                  annotationColName = "label")
#' 
signatureROCplot <- function(inputData, annotationData, signatureColNames,
                             annotationColName, scale = FALSE,
                             choose_colors = c("cornflowerblue", "grey24"),
                             name = "Signatures", nrow = NULL, ncol = NULL) {
  
  # Error catches and variable creation
  if (methods::is(inputData, "SummarizedExperiment")) {
    if (any(duplicated(signatureColNames))) {
      signatureColNames <- unique(signatureColNames)
    }
    if (!all(signatureColNames %in% colnames(SummarizedExperiment::colData(
      inputData)))) {
      stop("Signature column name not found in inputData.")
    }
    if (!all(annotationColName %in% colnames(SummarizedExperiment::colData(
      inputData)))) {
      stop("Annotation column name not found in inputData.")
    }
    annotationData <- data.frame(SummarizedExperiment::colData(
      inputData)[, annotationColName, drop = FALSE])
    inputData <-  data.frame(SummarizedExperiment::colData(
      inputData)[, signatureColNames, drop = FALSE])
  } else {
    if (ncol(annotationData) != 1) {
      stop("annotationData must have only one column.")
    }
    annotationColName <- colnames(annotationData)
  }
  if (length(annotationColName) != 1) {
    stop("You must specify a single annotation column name with which
    to color boxplots.")
  }
  if (!is.factor(annotationData[, 1])) {
    annotationData[, 1] <- as.factor(annotationData[, 1])
  }
  # The number of rows of annotation data should equal the
  # number of rows of the input data.
  if (nrow(annotationData) == nrow(inputData)) {
    if (!all(rownames(annotationData) == rownames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
  } else if (nrow(annotationData) == ncol(inputData)) {
    if (!all(rownames(annotationData) == colnames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
    inputData <- t(inputData)
  } else {
    stop("Annotation data and signature data does not match.")
  }
  
  # Begin plot creation
  pathwaydata <- t(inputData)
  if (scale) {
    pathwaydata <- t(scale(t(pathwaydata)))
  }
  inputdf <- data.frame(t(pathwaydata),
                        Group = annotationData[, 1])
  
  # If there are multiple annotationColNames, we will produce
  # an array of AUC plots as output.
  plot_dat <- NULL
  for (k in signatureColNames) {
    pred <- ROCit::rocit(inputdf[, k], inputdf$Group)
    plot_n <- as.data.frame(cbind(FPR = round(pred$FPR, 4), 
                                  TPR = round(pred$TPR, 4), Signature = k))
    plot_dat <- rbind(plot_dat, plot_n)
  }
  plot_dat$Signature <- paste(plot_dat$Signature)
  plot_dat$FPR <- as.numeric(paste(plot_dat$FPR))
  plot_dat$TPR <- as.numeric(paste(plot_dat$TPR))
  
  theplot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(ggplot2::aes(x = FPR, y = TPR, col = paste(
      "Empirical ROC curve")), size = 1) + 
    ggplot2::facet_wrap(~Signature, scales = 'free', 
                        nrow = nrow, ncol = ncol) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1, 
                                      col = choose_colors[2]), size = 1, 
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",
                  title = name) +
    ggplot2::scale_color_manual("", 
                                labels = c("Empirical ROC curve", "Chance line"),
                                values = c(choose_colors[1], choose_colors[2])) +
    ggplot2::theme(legend.position = "right")
  
  return(theplot)
}

#' Create an array of ROC plots with Confidence Interval Bands to Compare Signatures.
#' 
#' @inheritParams signatureROCplot
#' @param choose_colors a vector of length 3 definining the colors to be used
#' in the ROC plots. The default is \code{c("cornflowerblue", 
#' "grey50", "grey79")}.
#' @param name a character string giving the title of the boxplot. If 
#' \code{NULL}, the plot title will be 
#' \code{"ROC plots for Gene Signatures, <ci.lev>\% Confidence"}. 
#' The default is \code{NULL}.
#' @param ci.lev a number between 0 and 1 giving the desired level of 
#' confidence for computing ROC curve estimations.
#' 
#' @return An array of ROC plots.
#' 
#' @export
#' 
#' @examples 
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[-c(18, 22)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",  
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs)
#'                                  
#' # Create ROC plots with cocnfidence intervals
#' signatureROCplot_CI(prof_indian, signatureColNames = names(choose_sigs),
#'                     annotationColName = "label")
#' 
signatureROCplot_CI <- function(inputData, annotationData, signatureColNames,
                                annotationColName, scale = FALSE,
                                choose_colors = c("cornflowerblue", 
                                                  "grey50", "grey79"),
                                name = NULL, nrow = NULL, ncol = NULL,
                                ci.lev = 0.95) {
  
  # Error catches and variable creation
  if (methods::is(inputData, "SummarizedExperiment")) {
    if (any(duplicated(signatureColNames))) {
      signatureColNames <- unique(signatureColNames)
    }
    if (!all(signatureColNames %in% colnames(SummarizedExperiment::colData(
      inputData)))) {
      stop("Signature column name not found in inputData.")
    }
    if (!all(annotationColName %in% colnames(SummarizedExperiment::colData(
      inputData)))) {
      stop("Annotation column name not found in inputData.")
    }
    annotationData <- data.frame(SummarizedExperiment::colData(
      inputData)[, annotationColName, drop = FALSE])
    inputData <-  data.frame(SummarizedExperiment::colData(
      inputData)[, signatureColNames, drop = FALSE])
  } else {
    if (ncol(annotationData) != 1) {
      stop("annotationData must have only one column.")
    }
    annotationColName <- colnames(annotationData)
  }
  if (length(annotationColName) != 1) {
    stop("You must specify a single annotation column name with which
    to color boxplots.")
  }
  if (!is.factor(annotationData[, 1])) {
    annotationData[, 1] <- as.factor(annotationData[, 1])
  }
  # The number of rows of annotation data should equal the
  # number of rows of the input data.
  if (nrow(annotationData) == nrow(inputData)) {
    if (!all(rownames(annotationData) == rownames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
  } else if (nrow(annotationData) == ncol(inputData)) {
    if (!all(rownames(annotationData) == colnames(inputData))) {
      stop("Annotation data and signature data does not match.")
    }
    inputData <- t(inputData)
  } else {
    stop("Annotation data and signature data does not match.")
  }
  
  # Begin plot creation
  pathwaydata <- t(inputData)
  if (scale) {
    pathwaydata <- t(scale(t(pathwaydata)))
  }
  inputdf <- data.frame(t(pathwaydata),
                        Group = annotationData[, 1])
  
  # If there are multiple annotationColNames, we will produce
  # an array of AUC plots as output.
  
  plot_dat <- NULL
  for (k in signatureColNames) {
    pred <- ROCit::rocit(inputdf[, k], inputdf$Group)
    this.conf <- suppressWarnings(ROCit::ciROC(pred, level = ci.lev))
    if(any(is.nan(this.conf$UpperTPR))) {
      this.conf$UpperTPR[which(is.nan(this.conf$UpperTPR))] <- 1
    }
    
    plot_n <- as.data.frame(cbind(FPR = round(this.conf$FPR, 4), 
                                  TPR = round(this.conf$TPR, 4), 
                                  LowerTPR = round(this.conf$LowerTPR, 4),
                                  UpperTPR = round(this.conf$UpperTPR, 4),
                                  Signature = k))
    plot_dat <- rbind(plot_dat, plot_n)
  }
  
  plot_dat <- as.data.frame(apply(plot_dat, 2, paste))
  plot_dat[,1:4] <- as.data.frame(apply(plot_dat[,1:4], 2, as.numeric))
  
  if(is.null(name)) name <- paste("ROC plots for Gene Signatures, ", ci.lev, 
                                  "% Confidence", sep = "")
  
  theplot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LowerTPR, ymax = UpperTPR, 
                                      fill = "Confidence Interval"), 
                         show.legend = FALSE) +
    ggplot2::geom_line(ggplot2::aes(x = FPR, y = TPR, col = 
                                      "Empirical ROC curve"), size = 1) + 
    ggplot2::facet_wrap(~Signature, scales = 'free', 
                        nrow = nrow, ncol = ncol) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1, 
                                      col = choose_colors[2]), size = 1, 
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::scale_color_manual("", 
                                labels = c("Empirical ROC curve", "Chance line"),
                                values = c(choose_colors[1], choose_colors[2])) +
    ggplot2::scale_fill_manual("", 
                               labels = c("Confidence Interval"),
                               values = c(choose_colors[3])) +
    ggplot2::labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",
                  title = name) +
    ggplot2::theme(legend.position = "right")
  
  return(theplot)
}

#' Add SummarizedExperiment Assays to the Data Structure.
#' 
#' Given an input of a Summarized Experiment with a counts or CPM assay, This 
#' function creates additional assays for a gene expression count dataset 
#' to be used in further analysis.
#' 
#' @param SE_obj a \code{SummarizedExperiment} object containing gene expression
#' data. Required.
#' @param input_name a character string specifying the name of the assay to
#' be referenced for creating additional assays. Default is \code{"counts"}.
#' @param output_name a character string to concatenate to "log" when computing
#' a log assay. If \code{NULL}, then \code{input_name} will be substituted.
#' Only used if \code{log = TRUE}. Default is \code{NULL}.
#' @param log logical. Indicate whether an assay returned should be the log
#' of whichever assay is specified in \code{"output_name"}. If 
#' \code{counts_to_CPM = TRUE} as well, then a log CPM assay will also
#' be created. Default is \code{FALSE}. 
#' @param counts_to_CPM logical. This argument only applies if the 
#' \code{input_type} is a counts assay. If \code{TRUE}, then the output assays
#' will include a normalized CPM assay. If \code{counts_to_CPM = TRUE} as well, 
#' then a log CPM assay will also be created. Default is \code{TRUE}.  
#' @param  prior.counts an integer specifying the average count to be added to 
#' each observation to avoid taking the log of zero. Used only if 
#' \code{log = TRUE}. The default is \code{3}.
#' 
#' @return This function returns a \code{SummarizedExperiment} object with up 
#' to 3 additional assay types attached to the original inputted object.
#' \item{cpm}{Counts per million}
#' \item{logcpm}{Log counts per million}
#' \item{log_<output_name>}{Log of original inputted assay. 
#' \code{<output_name>} will be replaced }
#' 
#' 
#' @author Aubrey Odom
#' 
#' @export
#' 
#' @examples 
#' # Create a log assay of the original assay input
#' # TB_hiv dataset already has counts data
#' log_only <- mkAssay(TB_hiv, log = TRUE, counts_to_CPM = FALSE)
#' log_only
#' 
#' # Create a CPM assay
#' CPM_only <- mkAssay(TB_hiv)
#' CPM_only
#' 
#' # Create a logCPM, logcounts, and CPM assay
#' all_assays <- mkAssay(TB_hiv, log = TRUE)
#' all_assays
#' 

mkAssay <- function(SE_obj, input_name = "counts", output_name = NULL,
                    log = FALSE, counts_to_CPM = TRUE,
                    prior_counts = 3) {
  
  if (!(log | counts_to_CPM)) stop("At least counts_to_CPM or log must be TRUE.")
  
  # Identify the main assay to be referenced
  assay_main <- SummarizedExperiment::assay(SE_obj, input_name)
  if (is.null(output_name)) output_name <- input_name
  
  if (counts_to_CPM) {
    data.norm <- edgeR::DGEList(counts = assay_main)
    dge_data.norm <- edgeR::calcNormFactors(data.norm)
    CPM_data.norm <- edgeR::cpm(dge_data.norm, log = FALSE)
    SummarizedExperiment::assay(SE_obj, "cpm") <- CPM_data.norm
  } 
  
  if (log) {
    SummarizedExperiment::assay(
      SE_obj, paste("log", output_name, sep = "_")) <- log(assay_main + 
                                                             prior_counts)
  } 
  
  if (counts_to_CPM & log) {
    logCPM_data.norm <- edgeR::cpm(
      dge_data.norm, log = TRUE, prior.counts = prior_counts)
    SummarizedExperiment::assay(SE_obj, "log_cpm") <- logCPM_data.norm
  }
  
  return(SE_obj)
}















