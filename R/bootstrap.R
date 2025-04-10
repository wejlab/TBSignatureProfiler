globalVariables(c("BS_AUC", "FPR", "LowerTPR", "Signatures",
                  "TBsignatures", "TPR", "UpperTPR", "sigAnnotData"))

#' Bootstrap the AUC and conduct T-Tests for a collection of signatures.
#'
#' Run bootstrapping of the AUC and derive the p-value for a 2-sample t-test
#' for all signatures tested on a given dataset.
#'
#' @param SE_scored a \code{SummarizedExperiment} object with genes as the
#' row features and signature scores in the \code{colData}. There should also
#' be a column of annotation data. Required.
#' @param annotationColName a character string giving the column name in
#' \code{colData} that contains the annotation data. Required.
#' @param signatureColNames a vector of column names in the
#' \code{colData} that contain the signature score data. Required.
#' @param num.boot integer. The number of times to bootstrap the data. The
#' default is \code{100}.
#' @param pb.show logical for whether to show a progress bar while running code.
#' The default is \code{TRUE}.
#'
#' @return A list of length 5 returning a vector of p-values for a 2-sample
#' t-test, bootstrapped AUC values, an AUC value for using all scored values
#' for all signatures specified in \code{signatureColNames},
#' and values for the lower and upper bounds of a bootstrapped AUC confidence
#' interval using \code{pROC::roc()}.
#'
#' @export
#'
#' @examples
#'  # Run signature profiling
#'  choose_sigs <- list("madeupsig" = c("FCRL3", "OAS2", "IFITM3"))
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",
#'                                  algorithm = "ssGSEA",
#'                                  combineSigAndAlgorithm = TRUE,
#'                                  signatures = choose_sigs,
#'                                  parallel.sz = 1)
#'  # Bootstrapping
#'  booted <- bootstrapAUC(SE_scored = prof_indian, annotationColName = "label",
#'                         signatureColNames = names(choose_sigs), num.boot = 2)
#'  booted
#'
bootstrapAUC <- function(SE_scored, annotationColName, signatureColNames,
                         num.boot = 100, pb.show = TRUE) {
  annotationData <- SummarizedExperiment::colData(SE_scored)[annotationColName][, 1]
  if (!is.factor(annotationData)) annotationData <- as.factor(annotationData)

  # Create progress bar
  total <- length(signatureColNames)
  counter <- 0
  if (pb.show)  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

  # Initialize containers
  pvals <- aucs <- pROC_lower <- pROC_upper <- numeric(total)
  aucs_boot <- matrix(nrow = num.boot,
                      ncol = total,
                      dimnames = list(seq(num.boot),
                                      signatureColNames))

  for (i in signatureColNames) {
    which.sig <- which(signatureColNames == i)
    score <- SummarizedExperiment::colData(SE_scored)[i][, 1]
    # Conduct a 2-sample t-test on the scores and their
    # corresponding Tuberculosis group status
    pvals[which.sig] <- stats::t.test(score ~ annotationData)$p.value

    # Obtain AUC based on entire dataset
    auc <- ROCit::rocit(score, annotationData)$AUC
    aucs[which.sig] <- max(auc, 1 - auc)

    # Proceed with bootstrapping
    for (j in seq_len(num.boot)) {
      index <- sample(seq_along(score), replace = TRUE)
      tmp_score <- score[index]
      tmp_annotationData <- annotationData[index]
      if (length(unique(tmp_annotationData)) == 2) {
        pred <- ROCit::rocit(tmp_score, tmp_annotationData)
      } else {
        # If there's only one value to be predicted, resample
        index <- sample(seq_along(score), replace = TRUE)
        tmp_score <- score[index]
        tmp_annotationData <- annotationData[index]
        pred <- ROCit::rocit(tmp_score, tmp_annotationData)
      }
      tmp_auc <- max(pred$AUC, 1 - pred$AUC)
      aucs_boot[j, which.sig] <- tmp_auc
    }

    # pROC
    roc1 <- suppressMessages(pROC::roc(predictor = score,
                                       response = annotationData))
    conf <- pROC::ci.auc(roc1, method = "delong", progress = "none")
    pROC_lower[which.sig] <- round(conf[1], 4)
    pROC_upper[which.sig] <- round(conf[3], 4)

    # Update the progress bar
    counter <- counter + 1
    if (pb.show) utils::setTxtProgressBar(pb, counter)
  }

  if (pb.show) close(pb)
  return(list("P-values" = pvals, "Boot AUC Values" = aucs_boot,
              "Non-Boot AUC Values" = aucs,
              "pROC Lower" = pROC_lower, "pROC Upper" = pROC_upper))
}

#' Create a table of results for t-tests and bootstrapped AUCs for multiple scored signatures.
#'
#' This function collects the results of bootstrapping and t-tests for a scored
#' gene expression dataset and presents them using a JavaScript table with an
#' R interface, or as a \code{data.frame}.
#'
#' @inheritParams bootstrapAUC
#' @param output a character string indicating the table output format. Possible
#' values are \code{DataTable} and \code{data.frame}. The default is
#' \code{DataTable}.
#'
#' @export
#'
#' @return A JavaScript table with an R interface using the \code{DT} package.
#'
#' @examples
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[c(1, 2)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs,
#'                                  parallel.sz = 1,
#'                                  update_genes = FALSE)
#'  # Create table
#'  tableAUC(SE_scored = prof_indian, annotationColName = "label",
#'           signatureColNames = names(choose_sigs))
#'
#'  # Create data.frame object
#' h <-  tableAUC(SE_scored = prof_indian, annotationColName = "label",
#'                signatureColNames = names(choose_sigs),
#'                output = "data.frame",
#'                num.boot = 5)
#' head(h)
#'
tableAUC <- function(SE_scored, annotationColName, signatureColNames,
                     num.boot = 100, pb.show = TRUE, output = "DataTable") {
  BS.Results <- bootstrapAUC(SE_scored, annotationColName, signatureColNames,
                             num.boot, pb.show)
  pvals <- BS.Results[["P-values"]]
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  lowerAUC <- round(apply(aucs_boot, 2, stats::quantile,
                          probs = .025), 4)
  upperAUC <- round(apply(aucs_boot, 2, stats::quantile,
                          probs = .975), 4)
  return_table <- data.frame("Signature" = signatureColNames,
                             "P.value" = round(pvals, 4),
                             "neg10xLog(P.value)" =
                               round(-10 * log(pvals), 4),
                             "LowerAUC" = lowerAUC,
                             "AUC" = round(aucs, 4),
                             "UpperAUC" = upperAUC)
  return_table$Signature <- paste(return_table$Signature)
  rownames(return_table) <- c(seq(1, nrow(return_table)))
  
  # Create interactive table
  if (output == "DataTable") {
    return(DT::datatable(return_table[order(aucs, decreasing = TRUE), ],
                         options = list(scrollX = TRUE, pageLength = 10),
                         rownames = FALSE))
  } else if (output == "data.frame") {
    return(return_table)
  }
  
}

#' Create a comparison plot of boxplots for bootstrapped AUC values.
#'
#' Present the results of AUC bootstrapping for a collection of scored
#' signatures via boxplots.
#'
#' @inheritParams bootstrapAUC
#' @param num.boot an integer indicating the number of times to bootstrap the
#' data.
#' @param name a character string giving the overall title for the plot.
#' The default is \code{"Boxplot Comparison of Signature AUCs"}.
#' @param pb.show logical for whether to show a progress bar while running code.
#' Default is \code{TRUE}.
#' @param abline.col the color to be used for the dotted line at AUC = 0.5
#' (the chance line). The default is \code{"red"}.
#' @param fill.col the color to be used to fill the boxplots.
#' The default is \code{"white"}.
#' @param outline.col the color to be used for the boxplot outlines.
#' The default is \code{"black"}.
#' @param rotateLabels If \code{TRUE}, rotate labels. Default is \code{FALSE}.
#' @param violinPlot logical. Setting \code{violinPlot = TRUE} creates violin
#' plots in place of boxplots. The mean and +/- 1 standard deviation are added
#' to the violin plot interior for each signature.
#' The default is \code{FALSE}.
#'
#' @export
#'
#' @return A plot with side-by-side boxplots of bootstrapped AUC values for
#' each specified signature.
#'
#' @examples
#' # Run signature profiling
#'  choose_sigs <- TBsignatures[c("Zak_RISK_16", "Zhao_NANO_6")]
#'  prof_indian <- runTBsigProfiler(TB_indian[seq_len(25), ],
#'                                  useAssay = "logcounts",
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs,
#'                                  parallel.sz = 1)
#'  # Create boxplots
#'  compareBoxplots(prof_indian, annotationColName = "label",
#'                  signatureColNames = names(choose_sigs), rotateLabels = TRUE)
#'
compareBoxplots <- function(SE_scored, annotationColName, signatureColNames,
                            num.boot = 100,
                            name = "Boxplot Comparison of Signature AUCs",
                            pb.show = TRUE, abline.col = "red",
                            fill.col = "gray79", outline.col = "black",
                            rotateLabels = FALSE, violinPlot = FALSE) {
  # Obtain AUCs
  BS.Results <- bootstrapAUC(SE_scored, annotationColName, signatureColNames,
                             num.boot, pb.show)
  aucs_boot <- data.frame(BS.Results[["Boot AUC Values"]])
  colnames(aucs_boot) <- signatureColNames
  aucs <- apply(aucs_boot, 2, stats::median)

  # Create boxplots with ggplot2
  melted_data <- reshape2::melt(aucs_boot, measure.vars = signatureColNames,
                                variable.name = "Signatures",
                                value.name = "BS_AUC")
  melted_data$Signatures <- gdata::reorder.factor(
    x = melted_data$Signatures,
    new.order = names(sort(aucs)))
  melted_data <- melted_data[order(melted_data$Signatures), ]
  the_plot <- ggplot2::ggplot(data = melted_data, ggplot2::aes(Signatures,
                                                               BS_AUC))
  if (violinPlot) {
    the_plot <- the_plot +
      ggplot2::geom_violin(fill = fill.col, col = outline.col) +
      ggplot2::stat_summary(fun.data = function(x) {
        m <- mean(x)
        ymin <- m - stats::sd(x)
        ymax <- m + stats::sd(x)
        return(c(y = m, ymin = ymin, ymax = ymax))
      }, geom = "pointrange", color = outline.col)
  } else {
    the_plot <- the_plot +
      ggplot2::geom_boxplot(fill = fill.col, col = outline.col)
  }

  the_plot <- the_plot +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0.5, slope = 0,
                                      col = abline.col), linewidth = 1,
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(
                     margin = ggplot2::margin(r = 10))) +
    ggplot2::ggtitle(label = name) +
    ggplot2::ylab(label = "Bootstrapped AUCs")

  if (rotateLabels) the_plot <- the_plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1))

  return(the_plot)
}

#' Create an array of ROC plots to compare signatures.
#'
#' @inheritParams signatureBoxplot
#' @param choose_colors a \code{vector} of length 2 defining the colors to be used
#' in the ROC plots. The default is \code{c("cornflowerblue", "gray24")}.
#' @param signatureColNames a \code{vector} of the column names of \code{inputData}
#' that contain the signature data. If \code{inputData} is a
#' \code{SummarizedExperiment} object, these are the column names of the
#' object \code{colData}.
#'
#' @return An array of ROC plots.
#'
#' @export
#'
#' @examples
#' # Run signature profiling
#'  choose_sigs <- subset(TBsignatures,
#'                        !(names(TBsignatures) %in% c("Lee_4", "Roe_OD_4")))[c(1,2)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",
#'                                  algorithm = "ssGSEA",
#'                                  signatures = choose_sigs,
#'                                  parallel.sz = 1)
#'
#' # Create ROC plots
#' signatureROCplot(prof_indian, signatureColNames = names(choose_sigs),
#'                  annotationColName = "label")
#'
signatureROCplot <- function(inputData, annotationData,
                             signatureColNames,
                             annotationColName, scale = FALSE,
                             choose_colors = c("cornflowerblue", "gray24"),
                             name = "Signatures", nrow = NULL, ncol = NULL) {
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
    to create plots.")
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
    plot_n <- as.data.frame(cbind("FPR" = round(pred$FPR, 4),
                                  "TPR" = round(pred$TPR, 4), Signature = k))
    plot_dat <- rbind(plot_dat, plot_n)
  }
  plot_dat$Signature <- paste(plot_dat$Signature)
  plot_dat$FPR <- as.numeric(paste(plot_dat$FPR))
  plot_dat$TPR <- as.numeric(paste(plot_dat$TPR))

  theplot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(ggplot2::aes(x = FPR, y = TPR, col = paste(
      "Empirical ROC curve")), linewidth = 1) +
    ggplot2::facet_wrap(~Signature, scales = "free",
                        nrow = nrow, ncol = ncol) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1,
                                      col = choose_colors[2]), linewidth = 1,
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",
                  title = name) +
    ggplot2::scale_color_manual("",
                                labels = c("Empirical ROC curve", "Chance line"),
                                values = c(choose_colors[1], choose_colors[2])) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  return(theplot)
}

#' Create an array of ROC plots with confidence interval bands to compare signatures.
#'
#' @inheritParams signatureROCplot
#' @param choose_colors a vector of length 3 defining the colors to be used
#' in the ROC plots. The default is \code{c("cornflowerblue",
#' "gray50", "gray79")}.
#' @param signatureColNames a \code{vector} of the column names of \code{inputData}
#' that contain the signature data. If \code{inputData} is a
#' \code{SummarizedExperiment} object, these are the column names of the
#' object \code{colData}.
#' @param name a character string giving the title of the ROC plot. If
#' \code{NULL}, the plot title will be
#' \code{"ROC Plots for Gene Signatures, <ci.lev>\% Confidence"}.
#' The default is \code{NULL}.
#' @param ci.lev a number between 0 and 1 giving the desired level of
#' confidence for computing ROC curve estimations.
#' @param pb.show logical for whether to show a progress bar while running code.
#' The default is \code{TRUE}.
#'
#' @return An array of ROC plots.
#'
#' @export
#'
#' @examples
#' # Run signature profiling
#'
#'  choose_sigs <- TBsignatures[c(1, 2)]
#'  prof_indian <- runTBsigProfiler(TB_indian, useAssay = "logcounts",
#'                                  algorithm = "Zscore",
#'                                  signatures = choose_sigs,
#'                                  parallel.sz = 1)
#'
#' # Create ROC plots with confidence intervals
#' signatureROCplot_CI(prof_indian, signatureColNames = names(choose_sigs),
#'                     annotationColName = "label")
#'
signatureROCplot_CI <- function(inputData, annotationData, signatureColNames,
                                annotationColName, scale = FALSE,
                                choose_colors = c("cornflowerblue",
                                                  "gray50", "gray79"),
                                name = NULL, nrow = NULL, ncol = NULL,
                                ci.lev = 0.95, pb.show = TRUE) {

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
    to create plots.")
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

  # Create progress bar
  total <- length(signatureColNames)
  counter <- 0
  if (pb.show)  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

  for (k in signatureColNames) {
    pred <- ROCit::rocit(inputdf[, k], inputdf$Group)
    this.conf <- suppressWarnings(ROCit::ciROC(pred, level = ci.lev))
    if (any(is.nan(this.conf$UpperTPR))) {
      this.conf$UpperTPR[which(is.nan(this.conf$UpperTPR))] <- 1
    }

    plot_n <- as.data.frame(cbind("FPR" = round(this.conf$FPR, 4),
                                  "TPR" = round(this.conf$TPR, 4),
                                  "LowerTPR" = round(this.conf$LowerTPR, 4),
                                  "UpperTPR" = round(this.conf$UpperTPR, 4),
                                  "Signature" = k))
    plot_dat <- rbind(plot_dat, plot_n)

    # Update the progress bar
    counter <- counter + 1
    if (pb.show) utils::setTxtProgressBar(pb, counter)
  }

  plot_dat <- as.data.frame(apply(plot_dat, 2, paste))
  plot_dat[, seq(1, 4)] <- as.data.frame(apply(plot_dat[, seq(1, 4)], 2, as.numeric))

  if (is.null(name)) name <- paste("ROC Plots for Gene Signatures, ",
                                   ci.lev * 100,
                                  "% Confidence", sep = "")
  if (pb.show) close(pb)

  theplot <- ggplot2::ggplot(data = plot_dat, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LowerTPR, ymax = UpperTPR,
                                      fill = "Confidence Interval"),
                         show.legend = FALSE) +
    ggplot2::geom_line(ggplot2::aes(x = FPR, y = TPR, col =
                                      "Empirical ROC curve"), linewidth = 1) +
    ggplot2::facet_wrap(~Signature, scales = "free",
                        nrow = nrow, ncol = ncol) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0, slope = 1,
                                      col = choose_colors[2]), linewidth = 1,
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::scale_color_manual("",
                                labels = c("Empirical ROC curve", "Chance line"),
                                values = c(choose_colors[1], choose_colors[2])) +
    ggplot2::scale_fill_manual("",
                               labels = c("Confidence Interval"),
                               values = c(choose_colors[3])) +
    ggplot2::labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",
                  title = name) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  return(theplot)
}

#' Add SummarizedExperiment assays to the data structure.
#'
#' Given a \code{SummarizedExperiment} input with a counts or CPM assay, this
#' function creates additional assays for by computing the CPM, log, or both
#' of the input assay to be used in further analysis.
#'
#' @param SE_obj a \code{SummarizedExperiment} object containing count or CPM
#' data. Required.
#' @param input_name a character string specifying the name of the assay to
#' be referenced for creating additional assays. Default is \code{"counts"}.
#' @param output_name a character string to use in place of the \code{input_name}.
#' If \code{NULL}, then \code{input_name} will be substituted.
#' Default is \code{NULL}. See Return details for how names are altered.
#' @param log logical. Indicate whether an assay returned should be the log
#' of whichever assay is specified in \code{"output_name"}. If
#' \code{counts_to_CPM = TRUE} as well, then a log CPM assay will also
#' be created. Default is \code{FALSE}.
#' @param counts_to_CPM logical. This argument only applies if the
#' \code{input_type} is a counts assay. If \code{TRUE}, then the output assays
#' will include a normalized CPM assay. If \code{log = TRUE} as well,
#' then a log CPM assay will also be created. Default is \code{TRUE}.
#' @param  prior_counts a small integer specifying the average count to be added to
#' each observation to avoid taking the log of zero. Used only if
#' \code{log = TRUE}. The default is \code{3}.
#'
#' @return This function returns a \code{SummarizedExperiment} object with up
#' to 3 additional assay types attached to the original inputted object.
#' \item{\code{output_name}_cpm}{Counts per million}
#' \item{log_\code{output_name}_cpm}{Log counts per million}
#' \item{log_\code{output_name}}{Log of original input assay.}
#'
#' @author Aubrey Odom-Mabey
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
  if (!(log || counts_to_CPM)) {
    stop("At least counts_to_CPM or log must be TRUE.")
  } else if (!(input_name %in% names(SummarizedExperiment::assays(SE_obj)))) {
    stop("input_name must be an assay in the SE_obj")
  }
  # Identify the main assay to be referenced
  assay_main <- SummarizedExperiment::assay(SE_obj, input_name)
  if (is.null(output_name)) output_name <- input_name
  if (log) { # if no CPM, but yes log
    SummarizedExperiment::assay(
      SE_obj, paste("log", output_name, sep = "_")) <- log(assay_main + prior_counts)
  }
  if (counts_to_CPM) {
    data.norm <- edgeR::DGEList(counts = assay_main)
    dge_data.norm <- edgeR::calcNormFactors(data.norm)
    CPM_data.norm <- edgeR::cpm(dge_data.norm, log = FALSE)
    SummarizedExperiment::assay(
      SE_obj, paste(output_name, "cpm", sep = "_")) <- CPM_data.norm
    if (log) { # if log and counts_to_CPM
      # log of CPM
      logCPM_data.norm <- edgeR::cpm(
        dge_data.norm, log = TRUE, prior.counts = prior_counts)
      SummarizedExperiment::assay(
        SE_obj, paste("log", output_name, "cpm", sep = "_")) <- logCPM_data.norm
    }
  }
  return(SE_obj)
}
