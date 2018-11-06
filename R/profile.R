.myenv <- new.env(parent = emptyenv())

#' Run TB Signature Profiling
#'
#' Run profiling using a set of signatures and profiling algorithms
#'
#' @param input Input data as a SummarizedExperiment, data.frame, or matrix.
#' Required
#' @param useAssay When using SummarizedExperiment style input, the assay to
#' use for signature profiling. Required for SummarizedExperiment
#' @param signatures A list of signatures to run. If signatures is NULL, the
#' default set of signatures from the included TBsignatures data set is run.
#' For details run ?TBsignatures.
#' @param algorithm A list of algorithms to run. The default is GSVA, ssGSEA,
#' and ASSIGN.
#' @param combineSigAndAlgorithm Logical. If TRUE, the row names will be in
#' the form <signature>_<algorithm>. If FALSE, there will be a column named
#' 'algorithm' that lists which algorithm is used. If NULL, and one algorithm
#' was used, the algorithm will not be listed. The default is false.
#' @param assignDir If using ASSIGN, a directory to save intermediate ASSIGN
#' results. The default is NULL and intermediate results will not be saved.
#' @param outputFormat Output data format. Can be SummarizedExperiment, matrix,
#' or data.frame. The default is the same as input.
#' @param parallel.sz Number of processors to use when doing the calculations in
#' parallel for GSVA and ssGSEA. If 0, all cores are used. The default is 0.
#' @param ASSIGNiter The number of iterations in the ASSIGN MCMC. The default is
#' 100000.
#' @param ASSIGNburnin The number of burn-in iterations for ASSIGN. These
#' iterations are discarded when computing the posterior means of the model
#' parameters. The default is 50000.
#'
#' @return A data.frame of signature profiling results.
#'
#' @export
#'
#' @examples
#' #create some toy data to test ACS_COR_16 signature, with 5 samples with low
#' #and five samples with high expression of the signatures genes.
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$ACS_COR_16,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10))))
#' res <- runTBsigProfiler(mat_testdata, algorithm = "GSVA", parallel.sz = 1)
#' res["ACS_COR_16", ]
runTBsigProfiler <- function(input, useAssay = NULL,
                             signatures = NULL,
                             algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                             combineSigAndAlgorithm = FALSE, assignDir = NULL,
                             outputFormat = NULL, parallel.sz = 0,
                             ASSIGNiter = 100000, ASSIGNburnin = 50000) {
  if (is.null(signatures)){
    utils::data("TBsignatures", package = "TBSignatureProfiler", envir = .myenv)
    signatures <- .myenv$TBsignatures
  }
  runindata <- input
  if (methods::is(runindata, "SummarizedExperiment")){
    if (!is.null(useAssay)){
        runindata <- SummarizedExperiment::assay(input, useAssay)
    } else {
      stop("useAssay required for SummarizedExperiment Input")
    }
    if (!combineSigAndAlgorithm & length(algorithm) > 1){
      stop("SummarizedExperiment not supported with ",
           "combineSigAndAlgorithm FALSE.")
    }
  } else if (!is.null(useAssay)) {
    stop("useAssay only supported for SummarizedExperiment objects")
  }
  if (class(runindata) == "data.frame"){
    runindata <- as.matrix(runindata)
  }
  if (class(runindata) != "matrix"){
    stop("Invalid input data type. Accepted input formats are matrix, ",
         "data.frame, or SummarizedExperiment. Your input: ",
         as.character(class(input)))
  }
  if (!all(algorithm %in% c("GSVA", "ssGSEA", "ASSIGN"))){
    stop("Invalid algorithm. Supported algorithms are GSVA, ssGSEA, and ASSIGN")
  }

  gsvaRes <- NULL
  if ("GSVA" %in% algorithm){
    message("Running GSVA")
    gsvaRes <- GSVA::gsva(runindata, signatures,
                          parallel.sz = parallel.sz)
  }

  gsvaRes_ssgsea <- NULL
  if ("ssGSEA" %in% algorithm){
    message("Running ssGSEA")
    gsvaRes_ssgsea <- GSVA::gsva(runindata, signatures, method = "ssgsea",
                                 parallel.sz = parallel.sz)
  }

  assign_res <- NULL
  if ("ASSIGN" %in% algorithm){
    predir <- getwd()
    delete_intermediate <- FALSE
    if (is.null(assignDir)){
      assignDir <- tempfile("assign")
      dir.create(assignDir)
      delete_intermediate <- TRUE
    }
    setwd(assignDir)
    message("Running ASSIGN")
    for (i in names(signatures)){
      message(i)
      currlist <- signatures[i]
      currlist[[1]] <- currlist[[1]][currlist[[1]] %in% rownames(runindata)]
      if (length(currlist[[1]]) < 2){
        message("Not enough signature genes in ", i, ", so analysis will not run.")
      } else {
        if (!file.exists(i)){
          set.seed(1234)
          ASSIGN::assign.wrapper(testData = runindata, trainingLabel = NULL,
                                 geneList = currlist, adaptive_S = TRUE,
                                 iter = ASSIGNiter, burn_in = ASSIGNburnin,
                                 outputDir = i)
        } else {
          message("Result already exists. Delete to re-run.")
        }
      }
    }
    setwd(predir)
    assign_res <- as.matrix(t(ASSIGN::gather_assign_results(assignDir)))
    if (nrow(assign_res) == 0){
      assign_res <- NULL
    }
    if (delete_intermediate){
      unlink(assignDir, recursive = TRUE)
    } else {
      message("Intermediate ASSIGN results available at ", assignDir)
    }
  }

  sig_result <- NULL
  if (length(algorithm) == 1){
    if (!is.null(gsvaRes)){
      sig_result <- gsvaRes
    } else if (!is.null(gsvaRes_ssgsea)){
      sig_result <- gsvaRes_ssgsea
    } else if (!is.null(assign_res)){
      sig_result <- assign_res
    } else {
      stop("ERROR: all valid outputs are empty.")
    }
  } else {
    combined_res <- data.frame()
    #if combineSigAndAlgorithm TRUE, we can just concatenate the name to the sig
    #name and combine. Otherwise we need to add an 'algorithm' column and
    #SummarizedExperiment output is not supported.
    if (is.null(combineSigAndAlgorithm)){
      stop("You must choose whether or not to combine the ",
           "signature and algorithm name using combineSigAndAlgorithm.")
    } else if (combineSigAndAlgorithm){
      if (!is.null(gsvaRes)){
        rownames(gsvaRes) <- paste("GSVA", rownames(gsvaRes), sep = "_")
        combined_res <- gsvaRes
      }
      if (!is.null(gsvaRes_ssgsea)){
        rownames(gsvaRes_ssgsea) <- paste("ssGSEA", rownames(gsvaRes_ssgsea),
                                          sep = "_")
        if (nrow(combined_res) == 0){
          combined_res <- gsvaRes_ssgsea
        } else {
          combined_res <- rbind(combined_res, gsvaRes_ssgsea)
        }
      }
      if (!is.null(assign_res)){
        rownames(assign_res) <- paste("ASSIGN", rownames(assign_res), sep = "_")
        if (nrow(combined_res) == 0){
          combined_res <- assign_res
        } else {
          combined_res <- rbind(combined_res, assign_res)
        }
      }
    } else {
      if (!is.null(outputFormat)){
        if (outputFormat == "SummarizedExperiment") {
          stop("SummarizedExperiment not supported with",
               " combineSigAndAlgorithm FALSE.")
        }
      } else if (methods::is(input, "SummarizedExperiment")){
        stop("SummarizedExperiment not supported with ",
             "combineSigAndAlgorithm FALSE.")
      }
      if (!is.null(gsvaRes)){
        alg_col <- gsvaRes[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("GSVA", nrow(gsvaRes))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(gsvaRes)
        colnames(pathcol) <- "pathway"
        gsvaRes <- cbind(pathcol, alg_col, gsvaRes)
        combined_res <- gsvaRes
      }
      if (!is.null(gsvaRes_ssgsea)){
        alg_col <- gsvaRes_ssgsea[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("ssGSEA", nrow(gsvaRes_ssgsea))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(gsvaRes_ssgsea)
        colnames(pathcol) <- "pathway"
        gsvaRes_ssgsea <- cbind(pathcol, alg_col, gsvaRes_ssgsea)
        if (nrow(combined_res) == 0){
          combined_res <- gsvaRes_ssgsea
        } else {
          combined_res <- rbind(combined_res, gsvaRes_ssgsea)
        }
      }
      if (!is.null(assign_res)){
        alg_col <- assign_res[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("ASSIGN", nrow(assign_res))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(assign_res)
        colnames(pathcol) <- "pathway"
        assign_res <- cbind(pathcol, alg_col, assign_res)
        if (nrow(combined_res) == 0){
          combined_res <- assign_res
        } else {
          combined_res <- rbind(combined_res, assign_res)
        }
      }
      rownames(combined_res) <- NULL
    }
    sig_result <- as.matrix(combined_res)
  }

  if (is.null(outputFormat)){
    #output same as input
    if (class(input) %in% c("SummarizedExperiment", "SingleCellExperiment",
                             "SCtkExperiment")){
      SummarizedExperiment::colData(input) <-
        S4Vectors::cbind(SummarizedExperiment::colData(input),
                         S4Vectors::DataFrame(t(sig_result)))
      return(input)
    } else if (class(input) == "matrix") {
      return(sig_result)
    } else if (class(input) == "data.frame"){
      dfres <- data.frame(sig_result)
      colnames(dfres) <- colnames(sig_result)
      return(dfres)
    } else {
      stop("Output format error.")
    }
  } else if (outputFormat == "matrix"){
    return(sig_result)
  } else if (outputFormat == "data.frame"){
    dfres <- data.frame(sig_result)
    colnames(dfres) <- colnames(sig_result)
    return(dfres)
  } else if (outputFormat == "SummarizedExperiment"){
    outdata <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(data = runindata),
      colData = S4Vectors::DataFrame(t(sig_result)))
    return(outdata)
  } else {
    stop("Output format error.")
  }
}
