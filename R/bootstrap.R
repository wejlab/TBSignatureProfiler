#' Bootstrap the AUC and compute a p-value for each profiled signature
#'
#' Run bootstrapping of the AUC and derive the p-value for a t-test.
#' for all signatures tested on a given dataset.
#' @param SE_scored A SummarizedExperiment object with genes in the rows,
#' and signature scores in the columns. There should also be a column indicating TB status.
#' @param annotationData A factor of annotation data, indicating TB status.
#' @param SignatureColNames The column names in the SE_scored colData that contain 
#' the signature score data.
#' @param num.boot The number of times to bootstrap the data.
#' @param pb.show Logical for whether to show a progress bar while running code.
#' 
#' @return A list of length 3 returning a vector of p-values for a 2-sample t-test, bootstrapped AUC values,
#' and AUC value for using all scored values for all signatures specified in SignatureColNames.
#' 
#' @export
#' 
#' @examples 

bootstrapAUC <- function(SE_scored, annotationData, signatureColNames, num.boot = 100,
                         pb.show = TRUE){
  pvals <- aucs <- aucs_boot <- NULL
  
  if (!is.factor(annotationData)) annotationData <- as.factor(annotationData)
  
  # Create progress bar
  total <- length(names(signatureColNames))
  counter <- 0
  if(pb.show)  pb <- txtProgressBar(min = 0, max = total, 
                 label = "Computing bootstrapped AUC for each signature.",
                 style = 3)
  
  for (i in names(signatureColNames)){
    score <- SummarizedExperiment::colData(SE_scored)[i][,1]
    # Conduct a 2-sample t-test on the scores and their
    # corresponding Tuberculosis group status
    pvals <- c(pvals, t.test(score ~ annotationData)$p.value)
    
    # Obtain AUC based on entire dataset
    pred <- ROCit::rocit(score, annotationData)
    auc <- pred$AUC
    aucs <- c(aucs, max(auc, 1 - auc))
    
    # Proceed with bootstrapping
    tmp_aucs <- NULL
    for (j in 1:num.boot){
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

#' Create an interactive table of results for the bootstrapped AUC and p-values for a scored dataset.
#' 
#' Collect the results of bootstrapping and t-tests for a scored gene expression dataset
#' and present using a JavaScript table with an R interface.
#'
#' @param SE_scored A SummarizedExperiment object with genes in the rows,
#' and signature scores in the columns. There should also be a column indicating TB status.
#' @param annotationData A vector of annotation data, typically indicating TB status.
#' @param SignatureColNames The column names in the SE_scored colData that contain 
#' the signature score data.
#' @param num.boot The number of times to bootstrap the data.
#' @param pb_show Logical for whether to show a progress bar while running code.
#' 
#' @export
#' 
#' @return 
#' 
#' @example 

tableAUC <- function(SE_scored, annotationData, signatureColNames, num.boot, pb_show = TRUE){
  # Run the bootstrapping function
  BS.Results <- bootstrapAUC(SE_scored, annotationData, signatureColNames, num.boot, pb.show = pb_show)
  pvals <- BS.Results[["P-values"]]
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  
  # Create interactive table
  return(DT::datatable(cbind("Signatures" = names(signatureColNames), 
                             "P-values" = round(pvals, 4), 
                             "-10*Log(p-values)" = round(-10 * log(pvals), 4), 
                             "LowerAUC" = round(apply(aucs_boot, 2, quantile, probs = .05), 4), 
                             "AUCs" = round(aucs, 4), 
                             "UpperAUC"= round(apply(aucs_boot, 2, quantile, probs = .95), 4)
                             )[order(aucs, decreasing = T),], 
                       options = list(scrollX = T, pageLength = 30), rownames = F))
}

#' Create a comparison plot of boxplots for bootstrapped AUC values.
#' 
#' Present the results of AUC bootstrapping for the various scored signatures via boxplots.
#'
#' @param SE_scored A SummarizedExperiment object with genes in the rows,
#' and signature scores in the columns. There should also be a column indicating TB status.
#' @param annotationData A vector of annotation data, typically indicating TB status.
#' @param SignatureColNames The column names in the SE_scored colData that contain 
#' the signature score data.
#' @param num.boot The number of times to bootstrap the data.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex.
#' @param cex A numerical value giving the amount by which plotting text and symbols should be magnified 
#' relative to the default.
#' @param name A character string giving the overall title for the plot.
#' @param pb_show Logical for whether to show a progress bar while running code.
#' 
#' @export
#' 
#' @return A plot with side-by-side boxplots of bootstrapped AUC values for each specified signature.
#' 
#' @example 
compareBoxplots <- function(SE_scored, annotationData, signatureColNames, num.boot = 100,
                            cex.axis = 0.7, cex = 0.25,
                            name = "Boxplot Comparison of Signatures",
                            pb_show = TRUE){
  # Run the bootstrapping function
  BS.Results <- bootstrapAUC(SE_scored, annotationData, signatureColNames, num.boot, pb.show = pb_show)
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  
  colnames(aucs_boot) <- names(signatureColNames)
  boxplot(aucs_boot[,order(aucs)], las = 2, cex.axis = cex.axis, cex = cex,
          main = name, ylab = "Bootstrapped AUC Values")
  abline(h = 0.5, col = "red", lty = 2)
}

#' Create an array of ROC plots to compare signatures.
#' 
#' @param inputData Either a SummarizedExperiment object that contains the
#' signature data and annotation as colData columns, or a data.frame or matrix
#' of signature data. Required.
#' @param annotationData If inputData is a data.frame or matrix of signature
#' data, a data.frame or matrix of annotation data.
#' @param signatureColNames If inputData is a SummarizedExperiment, the column
#' names in colData that contain the signature data.
#' @param annotationColName If inputData is a SummarizedExperiment, the column
#' name in colData that contains the annotation data.
#' @param scale Scale the signature data. The default is FALSE.
#' 
#' @return An array of ROC plots.
#' 
#' @export
#' 
#' @example 
signatureROCplot <- function(inputData, annotationData, signatureColNames,
                             annotationColName, scale = FALSE) {
  
  # Error catches and variable creation
  if (methods::is(inputData, "SummarizedExperiment")){
    if (any(duplicated(signatureColNames))){
      signatureColNames <- unique(signatureColNames)
    }
    if (!all(signatureColNames %in% colnames(SummarizedExperiment::colData(inputData)))){
      stop("Signature column name not found in inputData.")
    }
    if (!all(annotationColName %in% colnames(SummarizedExperiment::colData(inputData)))){
      stop("Annotation column name not found in inputData.")
    }
    annotationData <- data.frame(SummarizedExperiment::colData(inputData)[, annotationColName, drop = FALSE])
    inputData <-  data.frame(SummarizedExperiment::colData(inputData)[, signatureColNames, drop = FALSE])
  } else {
    if (ncol(annotationData) != 1){
      stop("annotationData must have only one column.")
    }
    annotationColName <- colnames(annotationData)
  }
  if (length(annotationColName) != 1){
    stop("You must specify a single annotation column name to color boxplots by.")
  }
  if (!is.factor(annotationData[, 1])) {
    annotationData[, 1] <- as.factor(annotationData[, 1])
  }
  # The number of rows of annotation data should equal the
  # number of rows of the input data.
  if (nrow(annotationData) == nrow(inputData)){
    if (!all(rownames(annotationData) == rownames(inputData))){
      stop("Annotation data and signature data does not match.")
    }
  } else if (nrow(annotationData) == ncol(inputData)){
    if (!all(rownames(annotationData) == colnames(inputData))){
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

  par(mfrow = c(1,1))
  for (k in signatureColNames){
    pred <- ROCit::rocit(inputdf[, k], inputdf$Group)
    auc <- pred$AUC
    auc <- max(auc, 1 - auc)
    plot(pred, YIndex = F, values = F, 
         col = c("cornflowerblue", "grey29"))
    title(paste("ROC Curve,", k, "\n","(AUC = ", 
                round(auc, 3),")", sep = " ")) 
  }
  
}




