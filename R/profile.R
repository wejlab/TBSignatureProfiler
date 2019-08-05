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
#' profiling when \code{input} is a \code{SummarizedExperiment}. Required only for
#' input data of the class \code{SummarizedExperiment}. If null, the assay
#' used will be \code{"counts"}. The default is \code{NULL.}
#' @param signatures a \code{list} of signatures to run with their associated genes.
#' This list should be in the same format as \code{TBsignatures}, included in
#' the TBSignatureProfiler package. If \code{signatures = NULL}, the default set
#' of signatures \code{TBsignatures} list is used. For details, run
#' \code{?TBsignatures}. The default is \code{NULL}.
#' @param algorithm a vector of algorithms to run, or character string if only
#' one is desired. The default is \code{c("GSVA", "ssGSEA", "ASSIGN",
#' "PLAGE", "Zscore", "singscore")}.
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
#' Barbie,  D.A., Tamayo, P., Boehm, J.S., Kim, S.Y., Moody,	 S.E., Dunn, I.F., Schinzel, A.C.,
#' Sandy, P., Meylan, E., Scholl, C., et al. (2009).  Systematic RNA interference reveals
#' that oncogenic	KRAS-driven cancers require TBK1. \emph{Nature} \strong{462}, 108-112.
#' doi: \href{https://doi.org/10.1038/nature08460}{10.1038/nature08460}.
#'
#' Foroutan, M. et al. (2018). Single sample scoring of molecular phenotypes.
#' \emph{BMC Bioinformatics}, \bold{19}. doi:
#' \href{https://doi.org/10.1186/s12859-018-2435-4}{10.1186/s12859-018-2435-4}.
#'
#' Lee, E. et al. (2008). Inferring pathway activity toward precise disease
#' classification. \emph{PLoS Comp Biol}, 4(11):e1000217. doi:
#' \href{https://doi.org/10.1371/journal.pcbi.1000217}{10.1371/journal.pcbi.1000217}
#'
#' Shen, Y. et al. (2015). ASSIGN: context-specific genomic profiling of
#' multiple heterogeneous biological pathways. \emph{Bioinformatics}, \bold{31},
#' 1745-1753. doi:
#' \href{https://doi.org/10.1093/bioinformatics/btv031}{10.1093/bioinformatics/btv031}.
#'
#' Subramanian, A. et al. (2005). Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression profiles. \emph{PNAS},
#' \strong{102}, 15545-15550. doi:
#' \href{https://doi.org/10.1073/pnas.0506580102}{10.1073/pnas.0506580102}.
#'
#' Tomfohr, J. et al. (2005). Pathway level analysis of gene expression using
#' singular value decomposition. \emph{BMC Bioinformatics}, 6:225. doi:
#' \href{https://doi.org/10.1186/1471-2105-6-225}{10.1186/1471-2105-6-225}
#'
#' @source Profiling for the Z-Score, PLAGE, GSVA, ssGSEA algorithms are all
#' conducted with the Bioconductor \code{GSVA} package. Profiling for the
#' singscore algorithm is conducted with the Bioconductor \code{singscore}
#' package.
#'
#' @export
#'
#' @examples
#' ## Using a data.frame input/output
#'  # Create some toy data to test Zak_RISK_16 signature, using 5 samples with low
#'  # expression & five samples with high expression of the signatures genes.
#' df_testdata <- as.data.frame(rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$Zak_RISK_16,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10)))))
#' res <- runTBsigProfiler(input = df_testdata,
#'                         signatures = TBsignatures,
#'                         algorithm = c("GSVA", "ssGSEA"),
#'                         combineSigAndAlgorithm = FALSE,
#'                         parallel.sz = 1)
#' subset(res, res$pathway == "Zak_RISK_16")
#'
#' ## Using a SummarizedExperiment input/output
#'  # The TB_indian SummarizedExperiment data is included in the package.
#' GSVA_res <- runTBsigProfiler(input = TB_indian,
#'                              useAssay = "logcounts",
#'                              signatures = TBsignatures,
#'                              algorithm = c("GSVA"),
#'                              combineSigAndAlgorithm = FALSE,
#'                              parallel.sz = 1)
#' GSVA_res$Zak_RISK_16
#'
runTBsigProfiler <- function(input, useAssay = NULL,
                             signatures = NULL,
                             algorithm = c("GSVA", "ssGSEA", "ASSIGN",
                                           "PLAGE", "Zscore", "singscore"),
                             combineSigAndAlgorithm = FALSE, assignDir = NULL,
                             outputFormat = NULL, parallel.sz = 0,
                             ASSIGNiter = 100000, ASSIGNburnin = 50000) {
  if (is.null(signatures)) {
    # Override with global environment
    if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
      get("TBsignatures", envir = .GlobalEnv)
      signatures <- TBsignatures
    } else {
      utils::data("TBsignatures", package = "TBSignatureProfiler",
                  envir = .myenv)
      signatures <- .myenv$TBsignatures
    }
  }
  runindata <- input
  if (methods::is(runindata, "SummarizedExperiment")) {
    if (is.null(useAssay)){
      if ("counts" %in% names(SummarizedExperiment::assays(input))){
        useAssay <- "counts"
      } else {
        stop("useAssay required for SummarizedExperiment Input")
      }
    }
    runindata <- SummarizedExperiment::assay(input, useAssay)
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
  if (!all(algorithm %in% c("GSVA", "ssGSEA", "ASSIGN",
                            "PLAGE", "Zscore", "singscore"))) {
    stop("Invalid algorithm. Supported algorithms are
    GSVA, ssGSEA, PLAGE, Zscore, singscore, and ASSIGN")
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
  gsvaRes_PLAGE <- NULL
  if ("PLAGE" %in% algorithm) {
    message("Running PLAGE")
    gsvaRes_PLAGE <- GSVA::gsva(runindata, signatures, method = "plage",
                                parallel.sz = parallel.sz)
  }
  gsvaRes_Z <- NULL
  if ("Zscore" %in% algorithm) {
    message("Running Z-score profiling")
    gsvaRes_Z <- GSVA::gsva(runindata, signatures, method = "zscore",
                            parallel.sz = parallel.sz)
  }
  singscore_res <- NULL
  if ("singscore" %in% algorithm) {
    singscore_res <- matrix(ncol = ncol(runindata),
                            nrow = length(signatures),
                            dimnames = list(names(signatures),
                                            colnames(runindata)))
    rankDat <- singscore::rankGenes(runindata)
    for (sig in names(signatures)) {
      singscore_res[sig, ] <- suppressWarnings(singscore::simpleScore(
        rankData = rankDat,
        upSet = TBsignatures[[sig]],
        knownDirection = FALSE)$TotalScore)
    }
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
    } else if (!is.null(gsvaRes_Z)) {
      sig_result <- gsvaRes_Z
    } else if (!is.null(gsvaRes_PLAGE)) {
      sig_result <- gsvaRes_PLAGE
    } else if (!is.null(singscore_res)) {
      sig_result <- singscore_res
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
      if (!is.null(gsvaRes_PLAGE)) {
        rownames(gsvaRes_PLAGE) <- paste("PLAGE", rownames(gsvaRes_PLAGE),
                                         sep = "_")
        if (nrow(combined_res) == 0) {
          combined_res <- gsvaRes_PLAGE
        } else {
          combined_res <- rbind(combined_res, gsvaRes_PLAGE)
        }
      }
      if (!is.null(gsvaRes_Z)) {
        rownames(gsvaRes_Z) <- paste("Zscore", rownames(gsvaRes_Z),
                                     sep = "_")
        if (nrow(combined_res) == 0) {
          combined_res <- gsvaRes_Z
        } else {
          combined_res <- rbind(combined_res, gsvaRes_Z)
        }
      }
      if (!is.null(singscore_res)) {
        rownames(singscore_res) <- paste("singscore", rownames(singscore_res),
                                         sep = "_")
        if (nrow(combined_res) == 0) {
          combined_res <- singscore_res
        } else {
          combined_res <- rbind(combined_res, singscore_res)
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
      if (!is.null(gsvaRes_PLAGE)) {
        alg_col <- gsvaRes_PLAGE[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("PLAGE", nrow(gsvaRes_PLAGE))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(gsvaRes_PLAGE)
        colnames(pathcol) <- "pathway"
        gsvaRes_PLAGE <- cbind(pathcol, alg_col, gsvaRes_PLAGE)
        if (nrow(combined_res) == 0) {
          combined_res <- gsvaRes_PLAGE
        } else {
          combined_res <- rbind(combined_res, gsvaRes_PLAGE)
        }
      }
      if (!is.null(gsvaRes_Z)) {
        alg_col <- gsvaRes_Z[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("Z-Score", nrow(gsvaRes_Z))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(gsvaRes_Z)
        colnames(pathcol) <- "pathway"
        gsvaRes_Z <- cbind(pathcol, alg_col, gsvaRes_Z)
        if (nrow(combined_res) == 0) {
          combined_res <- gsvaRes_Z
        } else {
          combined_res <- rbind(combined_res, gsvaRes_Z)
        }
      }
      if (!is.null(singscore_res)) {
        alg_col <- singscore_res[, 1, drop = FALSE]
        colnames(alg_col) <- "algorithm"
        alg_col[, 1] <- rep("singscore", nrow(singscore_res))
        pathcol <- alg_col
        pathcol[, 1] <- rownames(singscore_res)
        colnames(pathcol) <- "pathway"
        singscore_res <- cbind(pathcol, alg_col, singscore_res)
        if (nrow(combined_res) == 0) {
          combined_res <- singscore_res
        } else {
          combined_res <- rbind(combined_res, singscore_res)
        }
      }
      rownames(combined_res) <- NULL
    }
    sig_result <- as.matrix(combined_res)
  }
  if (sum(names(signatures) %in% rownames(sig_result)) != length(signatures)) {
    absent <- subset(names(signatures), !(names(signatures) %in%
                                            rownames(sig_result)))
    if (length(algorithm) == 1) {
      warning(
        paste("No identifiers in the gene sets could be matched to the identifiers
          in the expression data for the following signatures: ",
              paste(absent, collapse = ", ")))
    }
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

#' Compare Scoring Algorithms on a Single Signature via Heatmap
#'
#' It may be useful to compare the results of scoring across several different
#' scoring algorithms via a method of visualization, such as a heatmap. The
#' \code{compareSigs} function allows the input of a data object and conducts
#' profiling on each signature desired, and outputting a heatmap for each
#' signature.
#'
#' @inheritParams runTBsigProfiler
#' @inheritParams signatureHeatmap
#' @param show.pb logical, whether warnings and other output
#' from the profiling should be suppressed (including progress bar output).
#' Default is \code{FALSE}.
#'
#' @return A heatmap for each signature specified comparing the enumerated
#' algorithms.
#'
#' @export
#'
#' @examples
#'
#' # Example using the TB_hiv data set, two signatures, and 3 algorithms
#' data("TB_hiv")
#' suppressWarnings({
#'   compareAlgs(TB_hiv, signatures = TBsignatures[c(1,2)],
#'               annotationColNames = "Disease",
#'               algorithm = c("GSVA", "ssGSEA", "PLAGE"),
#'               scale = TRUE, parallel.sz = 1)
#' })
#'
compareAlgs <- function (input, signatures = NULL, annotationColNames,
                         annotationData,
                         algorithm = c("GSVA", "ssGSEA", "ASSIGN", "PLAGE",
                                       "Zscore", "singscore"),
                         showColumnNames = TRUE,
                         showRowNames = TRUE, scale = FALSE,
                         useAssay = "counts",
                         colorSets = c("Set1", "Set2", "Set3", "Pastel1",
                                       "Pastel2", "Accent", "Dark2",
                                       "Paired"),
                         choose_color = c("blue", "gray95", "white"),
                         colList = list(),
                         show.pb = FALSE, parallel.sz = 0) {

  if (is.null(signatures)) {
    # Override with global environment
    if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
      get("TBsignatures", envir = .GlobalEnv)
      signatures <- TBsignatures
    } else {
      utils::data("TBsignatures", package = "TBSignatureProfiler",
                  envir = .myenv)
      signatures <- .myenv$TBsignatures
    }
  }
  for (sig in names(signatures)) {
    new.name <- paste("Scoring Methods for", sig)
    if (!show.pb) {
      utils::capture.output(scored <- runTBsigProfiler(input,
                                 useAssay = useAssay,
                                 combineSigAndAlgorithm = TRUE,
                                 signatures = signatures[sig],
                                 algorithm = algorithm,
                                 parallel.sz = parallel.sz))
    } else if (show.pb) {
      scored <- runTBsigProfiler(input, useAssay = useAssay,
                                 combineSigAndAlgorithm = TRUE,
                                 signatures = signatures[sig],
                                 algorithm = algorithm,
                                 parallel.sz = parallel.sz)
    }

    if (class(input) == "SummarizedExperiment") {
      already.there <- names(SummarizedExperiment::colData(input))
      col.names <- subset(names(SummarizedExperiment::colData(scored)),
                          !(names(SummarizedExperiment::colData(scored))
                            %in% already.there))
    } else col.names <- colnames(scored)

    return(signatureHeatmap(scored,
                            name = new.name,
                            annotationData = annotationData,
                            signatureColNames = col.names,
                            annotationColNames = annotationColNames,
                            scale = scale,
                            showColumnNames = showColumnNames,
                            showRowNames = showRowNames,
                            colorSets = colorSets, choose_color = choose_color,
                            colList = colList))
  }
}

