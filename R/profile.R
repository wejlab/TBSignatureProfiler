globalVariables(c("BS_AUC", "FPR", "LowerTPR", "Signatures",
                  "TBsignatures", "TPR", "UpperTPR", "sigAnnotData"))
.myenv <- new.env(parent = emptyenv())

# Helper functions for runTBsigProfiler

# Check whether a new sig list is in the working environment
check_sig_env <- function(signatures) {
  # Override with global environment
  if ("TBsignatures" %in% ls(envir = .GlobalEnv)) {
    get("TBsignatures", envir = .GlobalEnv)
    return(TBsignatures)
  } else {
    utils::data("TBsignatures", package = "TBSignatureProfiler",
                envir = .myenv)
    return(.myenv$TBsignatures)
  }
}

# Update gene names that are mogrified or outdated
update_genenames <- function(siglist) {
  requireNamespace("ggplot2", quietly = TRUE)
  if (!requireNamespace("HGNChelper", quietly = TRUE)) {
    stop("Package 'HGNChelper' required to update gene names.",
         "Either install or set update_genes = FALSE.")
  }
  newgenes <- suppressMessages(suppressWarnings(
    HGNChelper::checkGeneSymbols(siglist,
                                 unmapped.as.na = FALSE)))$Suggested.Symbol
  ind <- grep("//", newgenes)
  if (length(ind) != 0) newgenes[ind] <- strsplit(newgenes[ind],
                                                  " /// ")[[1]][1]
  # if(any(newgenes != siglist)) message("One or more gene names were altered.")
  return(newgenes)
}

# run scoring algorithms
score_algorithm <- function(runindata, signatures, alg, assignDir = NULL,
                            parallel.sz = 0,
                            ASSIGNiter = 100000, ASSIGNburnin = 50000,
                            ssgsea_norm = TRUE, combineSigAndAlgorithm,
                            one_alg) {
  if (alg %in% c("GSVA", "ssGSEA", "PLAGE", "Zscore")) {
    message("Running ", alg)
    algout <- GSVA::gsva(runindata, signatures, method = tolower(alg),
                         parallel.sz = parallel.sz, ssgsea.norm = ssgsea_norm)
  } else if (alg == "singscore") {
    message("Running singscore")
    algout <- matrix(ncol = ncol(runindata), nrow = length(signatures),
                     dimnames = list(names(signatures), colnames(runindata)))
    rankDat <- singscore::rankGenes(runindata)
    for (sig in names(signatures)) {
      algout[sig, ] <- suppressWarnings(singscore::simpleScore(
        rankData = rankDat,
        upSet = signatures[[sig]],
        knownDirection = FALSE)$TotalScore)
    }
  } else if (alg == "ASSIGN") {
    delete_intermediate <- FALSE
    if (is.null(assignDir)) {
      assignDir <- tempfile("assign")
      dir.create(assignDir)
      delete_intermediate <- TRUE
    } else if (!dir.exists(assignDir)) dir.create(assignDir)
    for (i in names(signatures)) {
      message("Running ASSIGN: ", i)
      currlist <- signatures[i]
      currlist[[1]] <- currlist[[1]][currlist[[1]] %in% rownames(runindata)]
      if (length(currlist[[1]]) < 2) {
        message("Not enough signature genes in ", i, "; cannot run.")
      } else {
        if (!file.exists(i)) {
          ASSIGN::assign.wrapper(testData = runindata, trainingLabel = NULL,
                                 geneList = currlist, adaptive_S = TRUE,
                                 iter = ASSIGNiter, burn_in = ASSIGNburnin,
                                 outputDir = file.path(assignDir, i))
        } else message("Result already exists. Delete to re-run.")
      }
    }
    algout <- as.matrix(t(ASSIGN::gather_assign_results(assignDir)))
    if (nrow(algout) == 0) algout <- NULL
    if (delete_intermediate) {
      unlink(assignDir, recursive = TRUE)
    } else {
      message("Intermediate ASSIGN results available at ", assignDir)
    }
  }
  if (one_alg) return(algout)
  if (length(combineSigAndAlgorithm) == 0) {
    stop("You must choose whether or not to combine the ",
         "signature and algorithm name using combineSigAndAlgorithm.")
  } else if (combineSigAndAlgorithm & !is.null(algout)) {
    scoremat <- algout
    rownames(scoremat) <- paste(alg, rownames(algout), sep = "_")
  } else if (!combineSigAndAlgorithm) {
    alg_col <- algout[, 1, drop = FALSE]
    colnames(alg_col) <- "algorithm"
    alg_col[, 1] <- rep(alg, nrow(algout))
    pathcol <- alg_col
    pathcol[, 1] <- rownames(algout)
    colnames(pathcol) <- "pathway"
    scoremat <- cbind(pathcol, alg_col, algout)
    rownames(scoremat) <- NULL
  }
  return(scoremat)
}

# Create proper output class for runTBsigProfiler
output_function <- function(input, signatures, algorithm, outputFormat,
                            runindata, sig_result) {
  if (is.null(outputFormat)) {
    # Output will be the same as input class
    if (class(input)[1] %in% c("SummarizedExperiment", "SingleCellExperiment",
                               "SCtkExperiment")) {
      SummarizedExperiment::colData(input) <-
        S4Vectors::cbind(SummarizedExperiment::colData(input),
                         S4Vectors::DataFrame(t(sig_result)))
      return(input)
    } else if (methods::is(input, "matrix")) {
      return(sig_result)
    } else if (methods::is(input, "data.frame")) {
      dfres <- data.frame(sig_result)
      colnames(dfres) <- colnames(sig_result)
      return(dfres)
    }
  } else if (outputFormat == "matrix") {
    return(sig_result)
  } else if (outputFormat == "data.frame") {
    dfres <- data.frame(sig_result)
    colnames(dfres) <- colnames(sig_result)
    return(dfres)
  } else if (outputFormat == "SummarizedExperiment") {
    attr(rownames(runindata), ".match.hash") <- NULL
    outdata <- SummarizedExperiment::SummarizedExperiment(
      assays = S4Vectors::SimpleList(data = runindata),
      colData = S4Vectors::DataFrame(t(sig_result)))
    return(outdata)
  } else stop("'OutputFormat' should be one of 'SummarizedExperiment', ",
              "'matrix', or 'data.frame'")
}


#' Run TB gene signature profiling.
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
#' of signatures \code{TBsignatures} list is used.
#' For details, run \code{?TBsignatures}.
#' If <2 genes in a signature are present in the sample, that signature will
#' not be evaluated and will not be present in the resulting SE object.
#' The default is \code{NULL}.
#' @param algorithm a vector of algorithms to run, or character string if only
#' one is desired. The default is \code{c("GSVA", "ssGSEA", "ASSIGN",
#' "PLAGE", "Zscore", "singscore")}.
#' @param combineSigAndAlgorithm logical, if \code{TRUE}, output row names will
#' be of the form <algorithm>_<signature>. It must be set to code{TRUE} if the
#' \code{ouputFormat} will be a SummarizedExperiment and
#' \code{length(algorithm) > 1}.
#' It will always be \code{FALSE} if only one algorithm is selected.
#' If \code{FALSE}, there will be a column named algorithm' that lists which
#' algorithm is used, and a column named 'pathway' that lists the signature
#' profiled. If \code{NULL}, and one algorithm was used, the algorithm will not
#' be listed. The default is \code{FALSE}.
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
#' @param ssgsea_norm logical, passed to \code{GSVA::gsva()}. When parameter
#' \code{algorithm = "ssgsea"},the profiler runs the SSGSEA method from
#' Barbie et al. (2009) normalizing the scores by the absolute difference
#' between the minimum and the maximum, as described in their paper.
#' When \code{ssgsea.norm = FALSE}, this last normalization step is skipped.
#' The default is \code{TRUE}.
#' @param update_genes logical, denotes whether gene names from \code{signatures}
#' and the rownames of \code{input} should be checked for accuracy using
#' \code{HGNChelper::checkGeneSymbols()}. The mapping assumes
#' genes are from humans and will keep unmapped genes as the original
#' input gene name. Default is \code{TRUE}.
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
#' Barbie, D.A., Tamayo, P., Boehm, J.S., Kim, S.Y., Moody, S.E., Dunn, I.F., Schinzel, A.C.,
#' Sandy, P., Meylan, E., Scholl, C., et al. (2009). Systematic RNA interference reveals
#' that oncogenic KRAS-driven cancers require TBK1. \emph{Nature} \strong{462}, 108-112.
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
#'                              paste0("sample", seq_len(10)))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", seq_len(100)),
#'                              paste0("sample", seq_len(10))))))
#' res <- runTBsigProfiler(input = df_testdata,
#'                         signatures = TBsignatures["Zak_RISK_16"],
#'                         algorithm = c("GSVA", "ssGSEA"),
#'                         combineSigAndAlgorithm = FALSE,
#'                         parallel.sz = 1)
#' subset(res, res$pathway == "Zak_RISK_16")
#'
#' ## Using a SummarizedExperiment input/output
#'  # The TB_indian SummarizedExperiment data is included in the package.
#' GSVA_res <- runTBsigProfiler(input = TB_indian,
#'                              useAssay = "logcounts",
#'                              signatures = TBsignatures["Zak_RISK_16"],
#'                              algorithm = c("GSVA"),
#'                              combineSigAndAlgorithm = FALSE,
#'                              parallel.sz = 1)
#' GSVA_res$Zak_RISK_16
runTBsigProfiler <- function(input, useAssay = NULL, signatures = NULL,
                             algorithm = c("GSVA", "ssGSEA", "ASSIGN",
                                           "PLAGE", "Zscore", "singscore"),
                             combineSigAndAlgorithm = FALSE, assignDir = NULL,
                             outputFormat = NULL, parallel.sz = 0,
                             ASSIGNiter = 100000, ASSIGNburnin = 50000,
                             ssgsea_norm = TRUE, update_genes = TRUE) {
  if (is.null(signatures)) signatures <- check_sig_env(signatures)
  acc_methods <- c("GSVA", "ssGSEA", "ASSIGN", "PLAGE", "Zscore", "singscore")
  if (!all(algorithm %in% acc_methods)) {
    stop("Invalid algorithm. Supported algorithms are: ",
         "GSVA, ssGSEA, PLAGE, Zscore, singscore, and ASSIGN")
  }
  runindata <- input
  if (methods::is(runindata, "SummarizedExperiment")) {
    if (is.null(useAssay)) {
      if ("counts" %in% names(SummarizedExperiment::assays(input))) {
        useAssay <- "counts"
      } else stop("useAssay required for SummarizedExperiment Input")
    }
    runindata <- SummarizedExperiment::assay(input, useAssay)
    if (!combineSigAndAlgorithm & length(algorithm) > 1) {
      stop("SummarizedExperiment not supported when ",
           "combineSigAndAlgorithm FALSE.")
    }
  } else if (!is.null(useAssay)) {
    stop("useAssay only supported for SummarizedExperiment objects")
  } else if (methods::is(runindata, "data.frame")) {
    runindata <- as.matrix(runindata)
  } else if (!methods::is(runindata, "matrix")) {
    stop("Invalid input data type. Accepted input formats are matrix, ",
         "data.frame, or SummarizedExperiment. Your input: ",
         as.character(class(input)))
  }
  if (update_genes) {
    message("Parameter update_genes is TRUE. Gene names will be updated.")
    signatures <- lapply(signatures, update_genenames)
    rownames(runindata) <- update_genenames(rownames(runindata))
  }
  # Remove signatures not present in the data
  sig_ind <- vapply(signatures,
                    function(x) sum(x %in% rownames(runindata)),
                    FUN.VALUE = numeric(1)) > 1
  if (any(!sig_ind)) message("The following signatures have <2 genes that",
                            " coincide with the genes in the given sample",
                            " and will not be scored: ",
                            paste0(names(signatures)[which(!sig_ind)],
                                   collapse = " "))
  signatures <- signatures[which(sig_ind)]
  if (length(algorithm) == 1) {
    sig_result <- score_algorithm(runindata, signatures, algorithm, assignDir,
                                  parallel.sz, ASSIGNiter, ASSIGNburnin,
                                  ssgsea_norm, combineSigAndAlgorithm,
                                  one_alg = TRUE)
    if (is.null(sig_result)) stop("All valid outputs are empty.")
  } else {
    combined <- lapply(algorithm, function(x) {
      score_algorithm(runindata, signatures, x, assignDir,
                      parallel.sz, ASSIGNiter, ASSIGNburnin,
                      ssgsea_norm, combineSigAndAlgorithm, one_alg = FALSE)
    })
    sig_result <- as.matrix(do.call(rbind, combined))
  }
  return(output_function(input, signatures, algorithm, outputFormat, runindata,
                         sig_result))
}

#' Compare scoring algorithms on a single signature via heatmap or boxplot.
#'
#' It may be useful to compare the results of scoring across several different
#' scoring algorithms via a method of visualization, such as a heatmap. The
#' \code{compareSigs} function allows the input of a SummarizedExperiment
#' data object and conducts
#' profiling on each signature desired, and outputting a heatmap or boxplot
#' for each signature.
#'
#' @inheritParams runTBsigProfiler
#' @inheritParams signatureHeatmap
#' @inheritParams compareBoxplots
#' @param input an input data object of the class \code{"SummarizedExperiment"}.
#' Required.
#' @param show.pb logical, whether warnings and other output
#' from the profiling should be suppressed (including progress bar output).
#' Default is \code{FALSE}.
#' @param output a character string specifying whether the outputted plot
#' should be a \code{"heatmap"} or \code{"boxplot"}. The default is
#' \code{"heatmap"}.
#'
#' @return A heatmap or boxplot for each signature specified comparing
#' the enumerated algorithms.
#'
#' @export
#'
#' @examples
#'
#' # Example using the TB_hiv data set, two signatures, and 3 algorithms
#' data("TB_hiv")
#' compareAlgs(TB_hiv, signatures = TBsignatures[c(1,2)],
#'             annotationColName = "Disease",
#'             algorithm = c("GSVA", "ssGSEA", "PLAGE"),
#'             scale = TRUE, parallel.sz = 1, output = "heatmap")
#'
compareAlgs <- function(input, signatures = NULL, annotationColName,
                        useAssay = "counts",
                        algorithm = c("GSVA", "ssGSEA", "ASSIGN", "PLAGE",
                                      "Zscore", "singscore"),
                        showColumnNames = TRUE,
                        showRowNames = TRUE, scale = FALSE,
                        colorSets = c("Set1", "Set2", "Set3", "Pastel1",
                                      "Pastel2", "Accent", "Dark2",
                                      "Paired"),
                        choose_color = c("blue", "gray95", "red"),
                        colList = list(),
                        show.pb = FALSE, parallel.sz = 0, output = "heatmap",
                        num.boot = 100, column_order = NULL) {
  if (output != "heatmap" & output != "boxplot") {
    stop("Output parameter must specify either 'heatmap' or 'boxplot'")
  }
  if (is.null(signatures)) signatures <- check_sig_env(signatures)
  for (sig in names(signatures)) {
    new.name <- paste("Scoring Methods for", sig)
    if (!show.pb) {
      utils::capture.output(scored <- runTBsigProfiler(
        input, useAssay = useAssay, combineSigAndAlgorithm = TRUE,
        signatures = signatures[sig], algorithm = algorithm,
        parallel.sz = parallel.sz))
    } else if (show.pb) {
      scored <- runTBsigProfiler(
        input, useAssay = useAssay, combineSigAndAlgorithm = TRUE,
        signatures = signatures[sig], algorithm = algorithm,
        parallel.sz = parallel.sz)
    }
    if (methods::is(input, "SummarizedExperiment")) {
      already.there <- names(SummarizedExperiment::colData(input))
      col.names <- subset(names(SummarizedExperiment::colData(scored)),
                          !(names(SummarizedExperiment::colData(scored))
                            %in% already.there))
    } else stop("Input must be a SummarizedExperiment object.")
    if (output == "heatmap") {
      return(signatureHeatmap(scored,
                              name = new.name,
                              signatureColNames = col.names,
                              annotationColNames = annotationColName,
                              scale = scale,
                              showColumnNames = showColumnNames,
                              showRowNames = showRowNames,
                              colorSets = colorSets, choose_color = choose_color,
                              colList = colList, split_heatmap = "none",
                              column_order = column_order))
    } else if (output == "boxplot") {
      return(compareBoxplots(scored,
                             annotationColName = annotationColName,
                             signatureColNames = col.names,
                             num.boot = num.boot, name = new.name,
                             pb.show = show.pb, rotateLabels = TRUE))
    }
  }
}
