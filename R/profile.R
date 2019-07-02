.myenv <- new.env(parent = emptyenv())

#' Run Tuberculosis Gene Signature profiling
#'
#' Using some subset of the signatures listed in \code{TBsignatures} and 
#' specified scoring algorithms, this function runs gene signature profiling 
#' on an input gene expression dataset. It allows for scores to be computed for 
#' these signatures which can be compared using various visualization tools also
#' provided in the TBSignatureProfiler package.
#'
#' @param input an input data object of the class \code{SummarizedExperiment}, 
#' \code{data.frame}, or \code{matrix} containing gene expression data. Required.
#' @param useAssay a character string specifying the assay to use for signature 
#' profiling when \code{input} is a SummarizedExperiment. Required only for input 
#' data of the class \code{SummarizedExperiment}.
#' @param signatures a \code{list} of signatures to run with their associated genes. 
#' This list should be in the same format as \code{TBsignatures}, included in 
#' the TBSignatureProfiler package. If \code{signatures = NULL}, the default set 
#' of signatures \code{TBsignatures} list is used. For details, run 
#' \code{?TBsignatures}.
#' @param algorithm a vector of algorithms to run, or character string if only 
#' one is desired. The default is \code{c("GSVA", "ssGSEA", "ASSIGN")}.
#' @param combineSigAndAlgorithm logical, not supported if \code{input} is a 
#' SummarizedExperiment object (in which case, the default is \code{TRUE}). 
#' For a matrix or data frame, if \code{TRUE}, the row names will be in the form 
#' <algorithm>_<signature>. If \code{FALSE}, there will be a column named
#' 'algorithm' that lists which algorithm is used, and a column named 'pathway'
#' that lists the signature profiled. If \code{NULL}, and one algorithm was used,
#' the algorithm will not be listed. The default is \code{FALSE}.
#' @param assignDir a character string naming a directory to save intermediate
#' ASSIGN results if \code{algorithm} specifies \code{"ASSIGN"}. The default is 
#' \code{NULL}, in which case intermediate results will not be saved.
#' @param outputFormat a character string specifying the output data format. 
#' Possible values are \code{"SummarizedExperiment"}, \code{"matrix"}, or 
#' \code{"data.frame"}. The default is to return the same type as the 
#' \code{input} object.
#' @param parallel.sz an integer identifying the number of processors to use 
#' when running the calculations in parallel for the GSVA and ssGSEA algorithms. 
#' If \code{parallel.sz = 0}, all cores are used. The default is \code{0}.
#' @param ASSIGNiter an integer indicating the number of iterations to use in 
#' the MCMC for the ASSIGN algorithm. The default is \code{100,000}.
#' @param ASSIGNburnin an integer indicating the number of burn-in iterations 
#' to use in the MCMC for the ASSIGN algorithm. These iterations are discarded 
#' when computing the posterior means of the model parameters. The default is 
#' \code{50,000}.
#'
#' @return A \code{SummarizedExperiment} object, \code{data.frame}, or 
#' \code{matrix} of signature profiling results. The returned object will be 
#' of the format specified in \code{outputFormat}. 
#' If \code{input} is a \code{SummarizedExperiment} and 
#' \code{outputFormat = "SummarizedExperiment"}, then the output will retain 
#' any input information stored in the input colData. In general, if 
#' \code{outputFormat = "SummarizedExperiment"} then columns in the \code{colData} 
#' will include the scores for each desired signature with samples on the rows. 
#' If \code{input} is a \code{data.frame} or \code{matrix}, then the returned 
#' object will have signatures on the rows and samples on the columns. 
#' 
#' @references 
#' Subramanian, A. et al. (2005). Gene set enrichment analysis: A knowledge-based 
#' approach for interpreting genome-wide expression profiles. \emph{PNAS},
#' \strong{102}, 15545-15550. doi:
#' \href{https://doi.org/10.1073/pnas.0506580102}{10.1073/pnas.0506580102}.
#' 
#' Hanzelmann, S., Castelo, R. and Guinney J. (2013). GSVA: gene set variation 
#' analysis for microarray and RNA-Seq data. \emph{BMC Bioinformatics}, 
#' \strong{14}. doi: 
#' \href{https://doi.org/10.1186/1471-2105-14-7}{10.1186/1471-2105-14-7}.
#' 
#' Shen, Y. et al. (2015). ASSIGN: context-specific genomic profiling of 
#' multiple heterogeneous biological pathways. \emph{Bioinformatics}, \bold{31},
#' 1745-1753. doi: 
#' \href{https://doi.org/10.1093/bioinformatics/btv031}{10.1093/bioinformatics/btv031}.
#'
#' @export
#'
#' @examples
#' ## Using a matrix input/output
#'  # Create some toy data to test ACS_COR_16 signature, using 5 samples with low 
#'  # expression & five samples with high expression of the signatures genes.
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$ACS_COR_16,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10))))
#' res <- runTBsigProfiler(input = mat_testdata,
#'                         signatures = TBsignatures, 
#'                         algorithm = c("GSVA", "ssGSEA"), 
#'                         combineSigAndAlgorithm = FALSE,
#'                         parallel.sz = 1)
#' res["ACS_COR_16", ]
#' 
#' ## Using a SummarizedExperiment input/output
#'  # The TB_indian SummarizedExperiment data is included in the package.
#' GSVA_res <- runTBsigProfiler(input = TB_indian,
#'                              useAssay = "logcounts",
#'                              signatures = TBsignatures,
#'                              algorithm = c("GSVA"),
#'                              combineSigAndAlgorithm = FALSE,
#'                              parallel.sz = 4)
#' GSVA_res$ACS_COR_16
#'  
runTBsigProfiler <- function(input, useAssay = NULL,
                             signatures = NULL,
                             algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                             combineSigAndAlgorithm = FALSE, assignDir = NULL,
                             outputFormat = NULL, parallel.sz = 0,
                             ASSIGNiter = 100000, ASSIGNburnin = 50000) {
  if (is.null(signatures)) {
    # Override with global environment
    if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
      get("TBsignatures", envir = .GlobalEnv)
      signatures <- TBsignatures
    } else {
      data("TBsignatures", package = "TBSignatureProfiler", envir = .myenv)
      signatures <- .myenv$TBsignatures
    }
  } 
  
  runindata <- input
  if (methods::is(runindata, "SummarizedExperiment")) {
    if (!is.null(useAssay)) {
        runindata <- SummarizedExperiment::assay(input, useAssay)
    } else {
      stop("useAssay required for SummarizedExperiment Input")
    }
    if (!combineSigAndAlgorithm & length(algorithm) > 1) {
      stop("SummarizedExperiment not supported with ",
           "combineSigAndAlgorithm FALSE.")
    }
  } else if (!is.null(useAssay)) {
    stop("useAssay only supported for SummarizedExperiment objects")
  }
  if (class(runindata) == "data.frame") {
    runindata <- as.matrix(runindata)
  }
  if (class(runindata) != "matrix") {
    stop("Invalid input data type. Accepted input formats are matrix, ",
         "data.frame, or SummarizedExperiment. Your input: ",
         as.character(class(input)))
  }
  if (!all(algorithm %in% c("GSVA", "ssGSEA", "ASSIGN"))) {
    stop("Invalid algorithm. Supported algorithms are GSVA, ssGSEA, and ASSIGN")
  }

  gsvaRes <- NULL
  if ("GSVA" %in% algorithm) {
    message("Running GSVA")
    gsvaRes <- GSVA::gsva(runindata, signatures,
                          parallel.sz = parallel.sz)
  }

  gsvaRes_ssgsea <- NULL
  if ("ssGSEA" %in% algorithm) {
    message("Running ssGSEA")
    gsvaRes_ssgsea <- GSVA::gsva(runindata, signatures, method = "ssgsea",
                                 parallel.sz = parallel.sz)
  }

  assign_res <- NULL
  if ("ASSIGN" %in% algorithm) {
    predir <- getwd()
    delete_intermediate <- FALSE
    if (is.null(assignDir)) {
      assignDir <- tempfile("assign")
      dir.create(assignDir)
      delete_intermediate <- TRUE
    }
    setwd(assignDir)
    message("Running ASSIGN")
    for (i in names(signatures)) {
      message(i)
      currlist <- signatures[i]
      currlist[[1]] <- currlist[[1]][currlist[[1]] %in% rownames(runindata)]
      if (length(currlist[[1]]) < 2) {
        message("Not enough signature genes in ", i, ", 
                so analysis will not run.")
      } else {
        if (!file.exists(i)) {
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
    if (delete_intermediate) {
      unlink(assignDir, recursive = TRUE)
    } else {
      message("Intermediate ASSIGN results available at ", assignDir)
    }
  }

  sig_result <- NULL
  if (length(algorithm) == 1) {
    if (!is.null(gsvaRes)) {
      sig_result <- gsvaRes
    } else if (!is.null(gsvaRes_ssgsea)) {
      sig_result <- gsvaRes_ssgsea
    } else if (!is.null(assign_res)) {
      sig_result <- assign_res
    } else {
      stop("ERROR: all valid outputs are empty.")
    }
  } else {
    combined_res <- data.frame()
    # if combineSigAndAlgorithm TRUE, we can just concatenate the name to the 
    # sig name and combine. Otherwise we need to add an 'algorithm' column and
    # SummarizedExperiment output is not supported.
    if (is.null(combineSigAndAlgorithm)) {
      stop("You must choose whether or not to combine the ",
           "signature and algorithm name using combineSigAndAlgorithm.")
    } else if (combineSigAndAlgorithm){
      if (!is.null(gsvaRes)){
        rownames(gsvaRes) <- paste("GSVA", rownames(gsvaRes), sep = "_")
        combined_res <- gsvaRes
      }
      if (!is.null(gsvaRes_ssgsea)) {
        rownames(gsvaRes_ssgsea) <- paste("ssGSEA", rownames(gsvaRes_ssgsea),
                                          sep = "_")
        if (nrow(combined_res) == 0) {
          combined_res <- gsvaRes_ssgsea
        } else {
          combined_res <- rbind(combined_res, gsvaRes_ssgsea)
        }
      }
      if (!is.null(assign_res)) {
        rownames(assign_res) <- paste("ASSIGN", rownames(assign_res), sep = "_")
        if (nrow(combined_res) == 0) {
          combined_res <- assign_res
        } else {
          combined_res <- rbind(combined_res, assign_res)
        }
      }
    } else {
      if (!is.null(outputFormat)) {
        if (outputFormat == "SummarizedExperiment") {
          stop("SummarizedExperiment not supported with",
               " combineSigAndAlgorithm FALSE.")
        }
      } else if (methods::is(input, "SummarizedExperiment")) {
        stop("SummarizedExperiment not supported with ",
             "combineSigAndAlgorithm FALSE.")
      }
      if (!is.null(gsvaRes)) {
        alg_col <- gsvaRes[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("GSVA", nrow(gsvaRes))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(gsvaRes)
        colnames(pathcol) <- "pathway"
        gsvaRes <- cbind(pathcol, alg_col, gsvaRes)
        combined_res <- gsvaRes
      }
      if (!is.null(gsvaRes_ssgsea)) {
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
      if (!is.null(assign_res)) {
        alg_col <- assign_res[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("ASSIGN", nrow(assign_res))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(assign_res)
        colnames(pathcol) <- "pathway"
        assign_res <- cbind(pathcol, alg_col, assign_res)
        if (nrow(combined_res) == 0) {
          combined_res <- assign_res
        } else {
          combined_res <- rbind(combined_res, assign_res)
        }
      }
      rownames(combined_res) <- NULL
    }
    sig_result <- as.matrix(combined_res)
  }

  if (sum(names(signatures) %in% rownames(sig_result)) != length(signatures)) {
    absent <- subset(names(signatures), !(names(signatures) %in% 
                                            rownames(sig_result)))
    warning(
      paste("No identifiers in the gene sets could be matched to the identifiers 
          in the expression data for the following signatures: ", 
            paste(absent, collapse = ", ")))
  }
  
  if (is.null(outputFormat)) {
    #output same as input
    if (class(input) %in% c("SummarizedExperiment", "SingleCellExperiment",
                             "SCtkExperiment")) {
      SummarizedExperiment::colData(input) <-
        S4Vectors::cbind(SummarizedExperiment::colData(input),
                         S4Vectors::DataFrame(t(sig_result)))
      return(input)
    } else if (class(input) == "matrix") {
      return(sig_result)
    } else if (class(input) == "data.frame") {
      dfres <- data.frame(sig_result)
      colnames(dfres) <- colnames(sig_result)
      return(dfres)
    } else {
      stop("Output format error.")
    }
  } else if (outputFormat == "matrix") {
    return(sig_result)
  } else if (outputFormat == "data.frame") {
    dfres <- data.frame(sig_result)
    colnames(dfres) <- colnames(sig_result)
    return(dfres)
  } else if (outputFormat == "SummarizedExperiment") {
    outdata <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(data = runindata),
      colData = S4Vectors::DataFrame(t(sig_result)))
    return(outdata)
  } else {
    stop("Output format error.")
  }
}


