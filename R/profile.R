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
#' #create some toy data to test Predict29 signature, with 5 samples with low
#' #and five samples with high expression of the signatures genes.
#' mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
#'                              dimnames = list(TBsignatures$Predict29,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10))))
#' res <- runTBsigProfiler(mat_testdata, algorithm = "GSVA", parallel.sz = 1)
#' res["Predict29", ]
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
      if (!dir.exists(assignDir)){
        dir.create(assignDir)
      } else {
        stop("Temp Directory Exists")
      }
      delete_intermediate <- TRUE
    }
    setwd(assignDir)
    message("Running ASSIGN")
    for (i in names(signatures)){
      message(i)
      currlist <- signatures[i]
      if (!all(signatures[[i]] %in% rownames(runindata))){
        message("Signature genes in ", i, " missing from input data:")
        message(paste(signatures[[i]][!(signatures[[i]] %in% rownames(runindata))], sep = "  ", collapse = ", "))
        message(i, " will be skipped")
      } else {
        if (!file.exists(paste0(i, "_ASSIGN"))){
          set.seed(1234)
          ASSIGN::assign.wrapper(testData = runindata, trainingLabel = NULL,
                                 geneList = currlist, adaptive_S = TRUE,
                                 iter = ASSIGNiter, burn_in = ASSIGNburnin,
                                 outputDir = paste0(i, "_ASSIGN"))
        }
      }
    }
    setwd(predir)
    assign_res <- as.matrix(t(ASSIGN::gather_assign_results(assignDir)))
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
    if (!is.null(gsvaRes)){
      return(gsvaRes)
    }
    if (!is.null(gsvaRes_ssgsea)){
      return(gsvaRes_ssgsea)
    }
    if (!is.null(assign_res)){
      return(assign_res)
    }
    sig_result <- combined_res
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
      return(data.frame(sig_result))
    } else {
      stop("Output format error.")
    }
  } else if (outputFormat == "matrix"){
    return(sig_result)
  } else if (outputFormat == "data.frame"){
    return(data.frame(sig_result))
  } else if (outputFormat == "SummarizedExperiment"){
    outdata <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(data = runindata),
      colData = S4Vectors::DataFrame(t(sig_result)))
    return(outdata)
  } else {
    stop("Output format error.")
  }
}
