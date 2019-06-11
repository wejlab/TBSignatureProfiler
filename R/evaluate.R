#' gene expression normalization, RLE normalization: relative log expression
#'
#' @param inputData a data.frame or matrix of gene expression count data.
#' Required.
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @return Normalized count data
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 10,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:10)))
#' res <- deseq2_norm_rle(inputTest)
deseq2_norm_rle <- function(inputData){
    scalingFac <- estimateSizeFactorsForMatrix(inputData)
    inputDataScaled <- inputData
    for (i in 1:ncol(inputData)){
        inputDataScaled[, i] <- inputData[, i] / scalingFac[i]
    }
    return(inputDataScaled)
}

#' LOOCV with logistic regression
#'
#' @param df a data.frame of gene expression count data.
#' Required.
#' @param targetVec binary vector indicating y
#' @import glmnet
#' @import caret
#' @importFrom ROCR prediction performance
#' @importFrom methods slot
#' @importFrom stats predict
#' @return AUC from LOOCV
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 10,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:10)))
#' inputTest <- as.data.frame(inputTest)
#' targetVec <- sample(c(0,1), replace=TRUE, size=10)
#' res <- LOOAUC_simple_multiple_noplot_one_df(inputTest, targetVec)

LOOAUC_simple_multiple_noplot_one_df <- function(df, targetVec){
  auc.vec <- c()
  nSample <- ncol(df)
  testPredictionClassVec <- c()
  testPredictionProbVec <- c()
  for (j in 1:nSample){
    train <- t(as.matrix(df[, -j]))
    test <- t(as.matrix(df[, j]))
    fit <- suppressWarnings(glmnet(train, targetVec[-j], family = "binomial"))
    testPredictionClassVec[j] <- suppressWarnings(predict(fit, type = "class", newx = test, s = 0))
    testPredictionProbVec[j] <- suppressWarnings(predict(fit, type = "response", newx = test, s = 0))
  }
  loo.pred <- suppressWarnings(prediction(testPredictionProbVec, targetVec))
  loo.perf <- suppressWarnings(performance(loo.pred, "tpr", "fpr"))
  auc <- suppressWarnings(performance(loo.pred, "auc"))
  auc <- unlist(slot(auc, "y.values"))
  aucRound <- round(auc, 4)
  auc.vec <- c(auc.vec, aucRound)
  # for other metric
  testPredictionClassVec <- as.numeric(testPredictionClassVec)
  cm <- suppressWarnings(confusionMatrix(as.factor(testPredictionClassVec), as.factor(targetVec)))
  output.list <- list()
  output.list[[1]] <- auc.vec
  output.list[[2]] <- cm$byClass
  output.list[[3]] <- testPredictionProbVec
  names(output.list) <- c("auc", "other", "prob")
  return(output.list)
}

#' Bootstrap LOOCV with logistic regression
#'
#' @param df a data.frame of gene expression count data.
#' Required.
#' @param target.vec binary vector indicating y
#' @param nboot number of bootstrap
#'
#' @return AUC from bootstrap LOOCV
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 20,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:20)))
#' inputTest <- as.data.frame(inputTest)
#' targetVec <- sample(c(0,1), replace=TRUE, size=20)
#' nboot <- 3
#' res <- Bootstrap_LOOCV_LR_AUC(inputTest, targetVec, nboot)
Bootstrap_LOOCV_LR_AUC <- function(df, target.vec, nboot){
  output.auc.vec <- c()
  output.other.df <- NULL
  for (i in 1:nboot){
    index.boot <- sample(1:ncol(df), ncol(df), replace = TRUE)
    df.tmp <- df[, index.boot]
    loo.output.list <- suppressWarnings(LOOAUC_simple_multiple_noplot_one_df(df.tmp, target.vec[index.boot]))
    output.auc.vec[i] <- loo.output.list[[1]]
    output.other.df <- rbind(output.other.df, loo.output.list[[2]])
  }

  output.list <- list()
  output.list[[1]] <- output.auc.vec
  output.list[[2]] <- as.data.frame(output.other.df)
  names(output.list) <- c("auc", "other")
  return(output.list)
}

#' Use logistic regression and Bootstrap LOOCV for signature evaluation.
#'
#' @param df.input a data.frame of gene expression count data.
#' Required.
#' @param target.vec.num numeric binary vector indicating y
#' @param signature.list a list of signatures
#' @param signature.name.vec vector of signature names
#' @param num.boot number of bootstrap
#' @importFrom gmodels ci
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics boxplot
#' @importFrom stats median
#' @return the AUC, sensitivity and specificity
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 20,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:20)))
#' inputTest <- as.data.frame(inputTest)
#' targetVec <- sample(c(0,1), replace=TRUE, size=20)
#' signature.list <- list(sig1 = c("gene1", "gene2", "gene3"),
#'                        sig2 = c("gene4", "gene5", "gene6"))
#' signature.name.vec <- c("sig1", "sig2")
#' num.boot <- 3
#' res <- SignatureQuantitative(inputTest,
#'                              targetVec,
#'                              signature.list,
#'                              signature.name.vec,
#'                              num.boot)
SignatureQuantitative <- function(df.input,
                                  target.vec.num,
                                  signature.list,
                                  signature.name.vec,
                                  num.boot = 100){

  df.list <- list()
  # create progress bar
  counter <- 0
  total <- length(signature.list)^2
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in 1:length(signature.list)){
    df.list[[i]] <- df.input[signature.list[[i]], ]
  }

  auc.result <- list()
  auc.result.ci <- list()
  sensitivity.ci <- list()
  specificity.ci <- list()
  
  for (i in 1:length(df.list)){
    boot.output.list <- suppressWarnings(Bootstrap_LOOCV_LR_AUC(df.list[[i]],
                                               target.vec.num,
                                               nboot = num.boot))
    # AUC
    auc.result[[i]] <- boot.output.list[[1]]
    result <- LOOAUC_simple_multiple_noplot_one_df(df.list[[i]], targetVec = target.vec.num)
    est <- result[[1]]
    ci.lower <- quantile(auc.result[[i]], probs = 0.05)
    ci.upper <- quantile(auc.result[[i]], probs = 0.95)
    st.error <- (1/(num.boot-1))*sum(auc.result[[i]] - mean(auc.result[[i]]))
    auc.result.ci[[i]] <- c("Estimate" = est, "CI lower" = ci.lower,
                            "CI upper" = ci.upper, "Std. Error" = st.error)
    names(auc.result)[i] <- signature.name.vec[i]
    names(auc.result.ci)[i] <- signature.name.vec[i]
    # sensitivity
    est2 <- result$other["Sensitivity"]
    ci.lower2 <- quantile(boot.output.list[[2]]$Sensitivity, probs = 0.05)
    ci.upper2 <- quantile(boot.output.list[[2]]$Sensitivity, probs = 0.95)
    st.error2 <- (1/(num.boot-1))*sum(boot.output.list[[2]]$Sensitivity - mean(boot.output.list[[2]]$Sensitivity)) 
    sensitivity.ci[[i]] <- c("Estimate" = est2, "CI lower" = ci.lower2,
                                   "CI upper" = ci.upper2, "Std. Error" = st.error2)
    names(sensitivity.ci)[i] <- signature.name.vec[i]
    # specificity
    est3 <- result$other["Specificity"]
    ci.lower3 <- quantile(boot.output.list[[2]]$Specificity, probs = 0.05)
    ci.upper3 <- quantile(boot.output.list[[2]]$Specificity, probs = 0.95)
    st.error3 <- (1/(num.boot-1))*sum(boot.output.list[[2]]$Specificity - mean(boot.output.list[[2]]$Specificity)) 
    specificity.ci[[i]] <- c("Estimate" = est3, "CI lower" = ci.lower3,
                                   "CI upper" = ci.upper3, "Std. Error" = st.error3)
    specificity.ci[[i]] <- suppressWarnings(ci(boot.output.list[[2]]$Specificity))
    names(specificity.ci)[i] <- signature.name.vec[i]
    
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }

  #boxplot
  m <- lapply(auc.result, median, na.rm = TRUE)
  o <- order(unlist(m))
  # pdf("boxplot.pdf", width = 13, height = 8)
  boxplot(auc.result[o],
          main = "Quantitative evaluation of signatures: AUC",
          cex.axis = 0.6,
          las = 2)
  # dev.off()

  # output df instead of list
  df.auc.ci <- data.frame(matrix(unlist(auc.result.ci),
                                 nrow = length(auc.result.ci),
                                 byrow = TRUE))
  colnames(df.auc.ci) <- names(auc.result.ci[[1]])
  rownames(df.auc.ci) <- signature.name.vec

  df.sensitivity.ci <- data.frame(matrix(unlist(sensitivity.ci),
                                         nrow = length(sensitivity.ci),
                                         byrow = TRUE))
  colnames(df.sensitivity.ci) <- names(sensitivity.ci[[1]])
  rownames(df.sensitivity.ci) <- signature.name.vec

  df.specificity.ci <- data.frame(matrix(unlist(specificity.ci),
                                         nrow = length(specificity.ci),
                                         byrow = TRUE))
  colnames(df.specificity.ci) <- names(specificity.ci[[1]])
  rownames(df.specificity.ci) <- signature.name.vec
  
  close(pb)

  return(list(df.auc.ci = df.auc.ci,
              df.sensitivity.ci = df.sensitivity.ci,
              df.specificity.ci = df.specificity.ci))
}

#' Add SummarizedExperiment assays to the Data Structure.
#' 
#' This function creates additional assays for a gene expression count dataset to be used in further analysis.
#' @param SE_obj SummarizedExperiment object containing gene expression count data. 
#' @param type vector or single character string indicating the type(s) of assay to add to SE_obj. 
#' If no value is specified, then CPM, log CPM, and log counts assays will all be created.
#' See details for further information on possible values.
#' @param  prior.counts integer specifying the average count to be added to each 
#' observation to avoid taking the log of zero. Used only if 'type' = "all" or "logCPM". 
#' If no value is input, default is 3.
#' 
#' @return This function returns a Summarized Experiment object with up 
#' to 3 additional assay types added to the original inputted object.
#' 
#' @details Here I will detail the specific objects... 
#' 
#' @author Aubrey Odom
#' 
#' @export
#' 
#' @examples 
#' ### You can add a log counts, log CPM, and CPM assay by default
#' testrun_all <- mkAssay(hivtb_data)
#' assays(testrun_all)
#' 
#' ### If a vector is passed to 'type' then all specified options will be added
#' testrun_CPM <- mkAssay(hivtb_data, type = c("CPM", "logCPM"))
#' assays(testrun_CPM)
#' 
#' ### Pass a single option to 'type'
#' testrun_lc <- mkAssay(hivtb_data, type = "logcounts")
#' assays(testrun_lc)

mkAssay <- function(SE_obj, type = "all", prior_counts = 3) {
  if (isFALSE(type %in% c("logcounts", "CPM", "logCPM", "all"))) {
    return(cat("\" The input for the assay type is invalid. Please specify one of the following options: 
                     \"logcounts\"
                     \"CPM\"
                     \"logCPM\"
                     \"all\"
                     "))
  } 
  
  if ("logcounts" %in% type) {
    SummarizedExperiment::assay(SE_obj, "logcounts") <- log(SummarizedExperiment::assays(SE_obj)$counts + 1)
  } 
  
  if ("CPM" %in% type) {
    data.norm <- edgeR::DGEList(counts = SummarizedExperiment::assay(SE_obj, "counts"))
    dge_data.norm <- edgeR::calcNormFactors(data.norm)
    CPM_data.norm <- edgeR::cpm(dge_data.norm, log = FALSE)
    SummarizedExperiment::assay(SE_obj, "cpm") <- CPM_data.norm
  } 
  
  if ("logCPM" %in% type) {
    data.norm <- edgeR::DGEList(counts = SummarizedExperiment::assay(SE_obj, "counts"))
    dge_data.norm <- edgeR::calcNormFactors(data.norm)
    logCPM_data.norm <- edgeR::cpm(dge_data.norm, log = TRUE, prior.counts = prior_counts)
    SummarizedExperiment::assay(SE_obj, "logcpm") <- logCPM_data.norm
  } 
  
  if ("all" %in% type) {
    # Create logcounts assay
    SummarizedExperiment::assay(SE_obj, "logcounts") <- log(SummarizedExperiment::assays(SE_obj)$counts + 1)
    
    # Create CPM assay
    data.norm <- edgeR::DGEList(counts = SummarizedExperiment::assay(SE_obj, "counts"))
    dge_data.norm <- edgeR::calcNormFactors(data.norm)
    CPM_data.norm <- edgeR::cpm(dge_data.norm, log = FALSE, prior.counts = prior_counts)
    SummarizedExperiment::assay(SE_obj, "cpm") <- CPM_data.norm
    
    # Create logCPM assay
    logCPM_data.norm <- edgeR::cpm(dge_data.norm, log = TRUE, prior.counts = prior_counts)
    SummarizedExperiment::assay(SE_obj, "logcpm") <- logCPM_data.norm
  }
  return(SE_obj)
}