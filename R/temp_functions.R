



## SE_scored is just the SE object produced from the profiler
## annotationData is a vector factor (?)
## signatureColNames is a list with elements that are signatures, 
## i.e. subset of the signature.list


# Depends on DT, ROCit, SummarizedExperiment




## BOOTSTRAP FUNCTION -----------------------------------------------

# Function that gives an interactive table of results
# for bootstrapped AUC and pvalues from t.test
bootstrapAUC <- function(SE_scored, annotationData, signatureColNames, num.boot = 100){
  pvals <- aucs <- aucs_boot <- NULL
  
  if (ncol(annotationData) != 1){
    stop("annotationData can only be one column.")
  } 
  if (!is.factor(annotationData)) annotationData <- as.factor(annotationData)
  
  # Create progress bar
  total <- length(names(signatureColNames))
  counter <- 0
  pb <- txtProgressBar(min = 0, max = total, 
                 label = "Computing bootstrapped AUC for each signature.",
                 style = 3)
  
  for (i in names(signatureColNames)){
    score <- SummarizedExperiment::colData(SE_scored)[i][,1]
    # Conduct a 2-sample t-test on the scores and their
    # corresponding Tuberculosis group status
    pvals <- c(pvals, t.test(score ~ annotationData)$p.value)
    
    # Obtain AUC based on entire dataset
    pred <- ROCit::rocit(SE_scored[, i], score)
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
    setTxtProgressBar(pb, counter)
  }
  
  close(pb)
  return(list("P-values" = pvals, "Boot AUC Values" = aucs_boot,
         "Non-Boot AUC Values" = aucs))
}


## TABLE FOR T-tests and AUC --------------------------------------

tableAUC <- function(SE_scored, annotationData, signatureColNames, num.boot){
  # Run the bootstrapping function
  BS.Results <- bootstrapAUC(SE_scored, annotationData, signatureColNames, num.boot)
  pvals <- BS.Results[["P-values"]]
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  
  # Create interactive table
  return(DT::datatable(cbind("Signatures" = names(signatureColNames), 
                             "P-values" = pvals, 
                             "-10*Log(p-values)" = round(-10 * log(pvals), 4), 
                             "LowerAUC" = round(apply(aucs_boot, 2, quantile, probs = .05), 4), 
                             "AUCs" = round(aucs, 4), 
                             "UpperAUC"= round(apply(aucs_boot, 2, quantile, probs = .95), 4)
                             )[order(aucs, decreasing = T),], 
                       options = list(scrollX = T, pageLength = 30), rownames = F))
}


## Boxplots for the bootstrapped AUC's ---------------------------

compareBoxplots <- function(SE_scored, annotationData, signatureColNames, num.boot = 100,
                            plot.col = 2, cex.axis = 0.7, cex = 0.25,
                            name = "Boxplot Comparison of Signatures"){
  # Run the bootstrapping function
  BS.Results <- bootstrapAUC(SE_scored, annotationData, signatureColNames, num.boot)
  aucs_boot <- BS.Results[["Boot AUC Values"]]
  aucs <- BS.Results[["Non-Boot AUC Values"]]
  
  colnames(aucs_boot) <- names(signatureColNames)
  boxplot(aucs_boot[,order(aucs)], las = 2, cex.axis = cex.axis, cex = cex,
          main = name, ylab = "Bootstrapped AUC Values")
  abline(h = 0.5, col = plot.col, lty = 2)
}

## AUC Plots ----------------------------------

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


