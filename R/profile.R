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
#'
#' @return A data.frame of signature profiling results.
#'
#' @export
#'
#' @examples
runTBsigProfiler <- function(input, useAssay = NULL,
                             signatures = NULL,
                             algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                             combineSigAndAlgorithm = FALSE, assignDir,
                             outputFormat = NULL) {
  if(is.null(signatures)){
    utils::data("TBsignatures", package = "TBSignatureProfiler", envir = .myenv)
    signatures <- .myenv$TBsignatures
  }
  runindata <- input
  if (!is.null(useAssay)){
    if (!(class(input) %in% c("SummarizedExperiment", "SingleCellExperiment",
                             "SCtkExperiment"))){
      stop("useAssay only supported for SummarizedExperiment objects")
    } else {
      runindata <- SummarizedExperiment::assay(input, useAssay)
    }
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
  if("GSVA" %in% algorithm){
    message("Running GSVA")
    gsvaRes <- GSVA::gsva(runindata,
                          signatures)
  }

  gsvaRes_ssgsea <- NULL
  if ("ssGSEA" %in% algorithm){
    message("Running ssGSEA")
    gsvaRes_ssgsea <- GSVA::gsva(runindata,
                                 signatures, method="ssgsea")
  }

  assign_res <- NULL
  if ("ASSIGN" %in% algorithm){
    message("Running ASSIGN")
  }

  sig_result <- NULL
  if(length(algorithm) == 1){
    if(!is.null(gsvaRes)){
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
    if(!is.null(gsvaRes)){
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

  if(is.null(outputFormat)){
    #output same as input
    if (class(input) %in% c("SummarizedExperiment", "SingleCellExperiment",
                             "SCtkExperiment")){
      SummarizedExperiment::colData(input) <-
        S4Vectors::cbind(SummarizedExperiment::colData(input),
                         S4Vectors::DataFrame(t(sig_result)))
      return(input)
    } else if (class(input) == "matrix") {
      return(sig_result)
    }
  } else if (outputFormat == "matrix"){
    return(sig_result)
  } else if(outputFormat == "data.frame"){
    return(data.frame(sig_result))
  } else if (outputFormat == "SummarizedExperiment"){
    outdata <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(data=runindata),
      colData=S4Vectors::DataFrame(t(sig_result)))
    return(outdata)
  }
  return(input)
}
