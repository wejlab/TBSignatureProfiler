#' Normalize Gene Expression Count Data.
#'
#' @param inputData a \code{data.frame} or \code{matrix} of gene expression
#' count data. Required.
#'
#' @return A \code{data.frame} or \code{matrix} of normalized count data.
#'
#' @export
#'
#' @examples
#' ## Example using the counts assay from a SummarizedExperiment
#' data_in <- SummarizedExperiment::assay(TB_indian, "counts")
#' res <- deseq2_norm_rle(data_in)
deseq2_norm_rle <- function(inputData){
    scalingFac <- DESeq2::estimateSizeFactorsForMatrix(inputData)
    inputDataScaled <- inputData
    for (i in 1:ncol(inputData)) {
        inputDataScaled[, i] <- inputData[, i] / scalingFac[i]
    }
    return(inputDataScaled)
}

#' Perform Leave-one-out CV with Logistic Regression.
#'
#' @param df a \code{data.frame} of gene expression count data. Required.
#' @param targetVec a binary vector of the response variable. Should be
#' the same number of rows as \code{df}. Required.
#'
#' @return A list of length 3 with elements
#' \item{auc}{The AUC from the LOOCV procedure.}
#' \item{byClass}{A vector containing the sensitivity, specificity, positive
#' predictive value, negative predictive value, precision, recall, F1,
#' prevalence, detection rate, detection prevalence and balanced accuracy.}
#' \item{prob}{A vector of the test prediction probabilities.}
#'
LOOAUC_simple_multiple_noplot_one_df <- function(df, targetVec){
  auc.vec <- c()
  nSample <- ncol(df)
  testPredictionClassVec <- c()
  testPredictionProbVec <- c()
  for (j in 1:nSample){
    train <- t(as.matrix(df[, -j]))
    test <- t(as.matrix(df[, j]))
    fit <- suppressWarnings(glmnet::glmnet(train, targetVec[-j],
                                           family = "binomial"))
    testPredictionClassVec[j] <- suppressWarnings(
      stats::predict(fit, type = "class", newx = test, s = 0))
    testPredictionProbVec[j] <- suppressWarnings(
      stats::predict(fit, type = "response", newx = test, s = 0))
  }
  loo.pred <- ROCit::rocit(testPredictionProbVec, targetVec)
  auc <- loo.pred$AUC
  aucRound <- round(auc, 4)
  auc.vec <- c(auc.vec, aucRound)
  testPredictionClassVec <- as.numeric(testPredictionClassVec)
  conf.mat <- suppressWarnings(
    caret::confusionMatrix(as.factor(testPredictionClassVec),
                                         as.factor(targetVec)))
  output.list <- list()
  output.list[[1]] <- auc.vec
  output.list[[2]] <- conf.mat$byClass
  output.list[[3]] <- testPredictionProbVec
  names(output.list) <- c("auc", "byClass", "prob")
  return(output.list)
}

#' Bootstrap on Leave-one-out CV with Logistic Regression.
#'
#' @param df a \code{data.frame} of gene expression count data. Required.
#' @param targetVec a binary vector of the response variable. Should be
#' the same number of rows as \code{df}. Required.
#' @param nboot an integer specifying the number of bootstrap iterations.
#'
#' @return A list of length 2 with elements \item{auc}{A vector the length of
#' \code{nboot} with the AUC from each bootstrap iteration.}
#' \item{byClass}{A dataframe with number of rows equal to \code{nboot}. Each
#' row contains the sensitivity, specificity, positive predictive
#' value, negative predictive value, precision, recall, F1, prevalence,
#' detection rate, detection prevalence and balanced accuracy for that
#' bootstrap iteration.}
#'
Bootstrap_LOOCV_LR_AUC <- function(df, targetVec, nboot){
  output.auc.vec <- c()
  output.byClass.df <- NULL
  for (i in 1:nboot){
    index.boot <- sample(1:ncol(df), ncol(df), replace = TRUE)
    df.tmp <- df[, index.boot]
    loo.output.list <- suppressWarnings(
      LOOAUC_simple_multiple_noplot_one_df(df.tmp, targetVec[index.boot]))
    output.auc.vec[i] <- loo.output.list[[1]]
    output.byClass.df <- rbind(output.byClass.df, loo.output.list[[2]])
  }

  output.list <- list()
  output.list[[1]] <- output.auc.vec
  output.list[[2]] <- as.data.frame(output.byClass.df)
  names(output.list) <- c("auc", "byClass")
  return(output.list)
}

#' Use Logistic Regression and Bootstrap LOOCV to Evaluate Signatures.
#'
#' This function takes as input a \code{data.frame} with genetic expression
#' count data, and uses a bootstrapped leave-one-out cross validation procedure
#' with logistic regression to allow for numeric and graphical comparison
#' across any number of genetic signatures.
#'
#' @inheritParams compareBoxplots
#' @param df.input a \code{data.frame} of gene expression count data. Required.
#' @param targetVec.num a numeric binary vector of the response variable.
#' The vector should be the same number of rows as \code{df}. Required.
#' @param signature.list a \code{list} of signatures to run with their
#' associated genes. This list should be in the same format as \code{TBsignatures},
#' included in the TBSignatureProfiler package. If \code{signature.list = NULL},
#' the default set of signatures \code{TBsignatures} list is used. For details,
#' run \code{?TBsignatures}.
#' @param signature.name.vec A vector specifying the names of the signatures
#' to be compared. This should be the same length as \code{signature.list}.
#' If \code{signature.name.vec = NULL}, the default set of signatures
#' \code{TBsignatures} list is used.
#' @param num.boot an integer specifying the number of bootstrap iterations.
#' @param pb.show logical. If \code{TRUE} then a progress bar for the
#' bootstrapping procedure will be displayed as output. The default is
#' \code{TRUE}.
#' @param name a character string giving a name for the outputted boxplot of
#' bootstrapped AUCs. The default is \code{"Quantitative Evaluation of
#' Signatures via Bootstrapped AUCs"}.
#'
#' @return the AUC, sensitivity and specificity
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 20,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:20)))
#' inputTest <- as.data.frame(inputTest)
#' targetVec <- sample(c(0,1), replace = TRUE, size = 20)
#' signature.list <- list(sig1 = c("gene1", "gene2", "gene3"),
#'                        sig2 = c("gene4", "gene5", "gene6"))
#' signature.name.vec <- c("sig1", "sig2")
#' num.boot <- 20
#' SignatureQuantitative(inputTest, targetVec.num = targetVec,
#'                       signature.list = signature.list,
#'                       signature.name.vec = signature.name.vec,
#'                       num.boot = num.boot)
SignatureQuantitative <- function(df.input, targetVec.num, signature.list = NULL,
                                  signature.name.vec = NULL, num.boot = 100,
                                  pb.show = TRUE) {

  if ((is.null(signature.name.vec) & !is.null(signature.list)
       | (!is.null(signature.name.vec) & is.null(signature.list)))){
    stop("Please specify arguments for both signature.list and
         signature.name.vec, or leave them both empty to use
         TBsignatures as the list of signatures for profiling.")
  } else if (is.null(signature.list) & is.null(signature.name.vec)){
    if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
      get("TBsignatures", envir = .GlobalEnv)
    }
    signature.list <- TBsignatures
    signature.name.vec <- names(signature.list)
  }

  if (length(signature.list) != length(signature.name.vec)){
    stop("The inputs signature.list and signature.name.vec are not the same
         length.")
  }

  df.list <- list()
  # progress bar
  counter <- 0
  total <- length(signature.list)
  if (pb.show) pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

  for (i in 1:length(signature.list)) {
    df.list[[i]] <- df.input[signature.list[[i]], ]
  }

  auc.result <- list()
  auc.result.ci <- list()
  sensitivity.ci <- list()
  specificity.ci <- list()

  for (i in 1:length(df.list)) {
    boot.output.list <- suppressWarnings(Bootstrap_LOOCV_LR_AUC(df.list[[i]],
                                                                targetVec.num,
                                                                nboot = num.boot))
    # AUC
    auc.result[[i]] <- boot.output.list[[1]]
    result <- LOOAUC_simple_multiple_noplot_one_df(df.list[[i]],
                                                   targetVec = targetVec.num)
    est <- result[[1]]
    ci.lower <- stats::quantile(auc.result[[i]], probs = 0.05)
    ci.upper <- stats::quantile(auc.result[[i]], probs = 0.95)
    st.error <- (1 / (num.boot - 1)) * sum(auc.result[[i]] -
                                             mean(auc.result[[i]]))
    auc.result.ci[[i]] <- c("Estimate" = est, "CI lower" = ci.lower,
                            "CI upper" = ci.upper, "Std. Error" = st.error)
    names(auc.result)[i] <- signature.name.vec[i]
    names(auc.result.ci)[i] <- signature.name.vec[i]
    # sensitivity
    est2 <- result$byClass["Sensitivity"]
    ci.lower2 <- stats::quantile(boot.output.list[[2]]$Sensitivity,
                                 probs = 0.05)
    ci.upper2 <- stats::quantile(boot.output.list[[2]]$Sensitivity,
                          probs = 0.95)
    st.error2 <- (1 / (num.boot - 1)) *
      sum(boot.output.list[[2]]$Sensitivity -
            mean(boot.output.list[[2]]$Sensitivity))
    sensitivity.ci[[i]] <- c("Estimate" = est2, "CI lower" = ci.lower2,
                             "CI upper" = ci.upper2, "Std. Error" = st.error2)
    names(sensitivity.ci)[i] <- signature.name.vec[i]
    est3 <- result$byClass["Specificity"]
    ci.lower3 <- stats::quantile(boot.output.list[[2]]$Specificity,
                                 probs = 0.05)
    ci.upper3 <- stats::quantile(boot.output.list[[2]]$Specificity,
                          probs = 0.95)
    st.error3 <- (1 / (num.boot - 1)) *
      sum(boot.output.list[[2]]$Specificity -
            mean(boot.output.list[[2]]$Specificity))
    specificity.ci[[i]] <- c("Estimate" = est3, "CI lower" = ci.lower3,
                             "CI upper" = ci.upper3, "Std. Error" = st.error3)
    specificity.ci[[i]] <- suppressWarnings(
      gmodels::ci(boot.output.list[[2]]$Specificity))
    names(specificity.ci)[i] <- signature.name.vec[i]

    counter <- counter + 1
    if (pb.show) utils::setTxtProgressBar(pb, counter)
  }

  # output data.frame instead of list
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

  if (pb.show) close(pb)

  return(list(df.auc.ci = df.auc.ci,
              df.sensitivity.ci = df.sensitivity.ci,
              df.specificity.ci = df.specificity.ci))
}

#' Create a Boxplot Using Logistic Regression and Bootstrap LOOCV to Evaluate Signatures.
#'
#' This function takes as input a \code{data.frame} with genetic expression
#' count data, and uses a bootstrapped leave-one-out cross validation procedure
#' with logistic regression to allow for numeric and graphical comparison
#' across any number of genetic signatures. It creates a boxplot of bootstrapped
#' AUC values.
#'
#' @inheritParams signatureBoxplot
#' @inheritParams SignatureQuantitative
#' @inheritParams compareBoxplots
#' @param name a character string giving a name for the outputted boxplot of
#' bootstrapped AUCs. The default is \code{"Signature Evaluation:
#' Bootstrapped AUCs"}.
#'
#' @return a boxplot comparing the bootstrapped AUCs of inputted signatures
#'
#' @export
#'
#' @examples
#' inputTest <- matrix(rnorm(1000), 100, 20,
#'                     dimnames = list(paste0("gene", 1:100),
#'                                     paste0("sample", 1:20)))
#' inputTest <- as.data.frame(inputTest)
#' targetVec <- sample(c(0,1), replace = TRUE, size = 20)
#' signature.list <- list(sig1 = c("gene1", "gene2", "gene3"),
#'                        sig2 = c("gene4", "gene5", "gene6"))
#' signature.name.vec <- c("sig1", "sig2")
#' num.boot <- 20
#' plotQuantitative(inputTest, targetVec.num = targetVec,
#'                  signature.list = signature.list,
#'                  signature.name.vec = signature.name.vec,
#'                  num.boot = num.boot, rotateLabels = FALSE)

plotQuantitative <- function(df.input, targetVec.num, signature.list = NULL,
                             signature.name.vec = NULL, num.boot = 100,
                             pb.show = TRUE, name =
                               "Signature Evaluation: Bootstrapped AUCs",
                             fill.col = "white", outline.col = "black",
                             abline.col = "red", rotateLabels = FALSE) {
  if ((is.null(signature.name.vec) & !is.null(signature.list)
       | (!is.null(signature.name.vec) & is.null(signature.list)))){
    stop("Please specify arguments for both signature.list and
         signature.name.vec, or leave them both empty to use
         TBsignatures as the list of signatures for profiling.")
  } else if (is.null(signature.list) & is.null(signature.name.vec)){
    if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
      get("TBsignatures", envir = .GlobalEnv)
    }
    signature.list <- TBsignatures
    signature.name.vec <- names(signature.list)
  }

  if (length(signature.list) != length(signature.name.vec)){
    stop("The inputs signature.list and signature.name.vec are not the same
         length.")
  }

  df.list <- list()
  # progress bar
  counter <- 0
  total <- length(signature.list)
  if (pb.show) pb <- utils::txtProgressBar(min = 0, max = total, style = 3)

  for (i in 1:length(signature.list)) {
    df.list[[i]] <- df.input[signature.list[[i]], ]
  }

  auc.result <- list()

  for (i in 1:length(df.list)) {
    boot.output.list <- suppressWarnings(Bootstrap_LOOCV_LR_AUC(df.list[[i]],
                                                                targetVec.num,
                                                                nboot = num.boot))
    # AUC
    auc.result[[i]] <- boot.output.list[[1]]
    counter <- counter + 1
    if (pb.show) utils::setTxtProgressBar(pb, counter)
  }

  # Boxplot
  auc.result <- data.frame(matrix(unlist(auc.result),
                                  ncol = length(signature.list),
                                  dimnames = list(c(), names(signature.list))))
  aucs <- apply(auc.result, 2, stats::median)

  melted_data <- reshape2::melt(auc.result,
                                measure.vars = names(signature.list),
                                variable.name = "Signatures",
                                value.name = "BS_AUC")
  melted_data$Signatures <- DescTools::reorder.factor(
    x = melted_data$Signatures,
    new.order = names(sort(aucs)))
  melted_data <- melted_data[order(melted_data$Signatures), ]
  the_plot <- ggplot2::ggplot(data = melted_data, ggplot2::aes(Signatures,
                                                               BS_AUC)) +
    ggplot2::geom_boxplot(fill = fill.col, col = outline.col) +
    ggplot2::geom_abline(ggplot2::aes(intercept = 0.5, slope = 0,
                                      col = abline.col), size = 1,
                         linetype = "dashed", show.legend = FALSE) +
    ggplot2::ggtitle(label = name) +
    ggplot2::ylab(label = "Bootstrapped AUCs") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(
                     margin = ggplot2::margin(r = 10)))

  if (rotateLabels) {
    the_plot <- the_plot + ggplot2::theme(axis.text.x = ggplot2::
                                          element_text(angle = 90, hjust = 1))
  }

  if (pb.show) close(pb)

  return(the_plot)
}
