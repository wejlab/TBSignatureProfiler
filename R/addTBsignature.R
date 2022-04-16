#' Introduce a new signature into the TBSignatureProfiler.
#'
#' This function allows users to integrate new signatures into the TBSP
#' with a function that updates the TBsignatures, TBcommon,
#' sigAnnotData and common_sigAnnotData objects. Users that wish to use this
#' function should do so with the downloaded package as a working directory,
#' and not as a casual package function. This function does not complete all
#' required updates to the package for a signature to be full added; users
#' should check the vignette "Submitting Signatures to the TBSP Package"
#' on the
#' [TBSP website](https://compbiomed.github.io/TBSignatureProfiler-docs/.)
#' for a walkthrough of this complete process.
#' Also note that this function only adds one signature at a time,
#' and must me run multiple times to add subsequent signatures.
#'
#' @param sigsymbols a \code{character vector} of the gene symbols that compose
#' the signature to be added. Required.
#' @param authname a \code{character} string containing the last name of the primary
#' author of the publication where the signature was first identified.
#' If spaces are present, omit them, and use proper capitalization. Required.
#' @param signame_common a \code{character} string of the alternate name of the
#' signature given by the publication, if it exists. If \code{NULL}, no
#' assigned name is assumed to exist. Default is \code{null}.
#' @param sigtype a \code{character} string that gives the context that the
#' signature was developed under. Most commonly,
#' it will distinguish TB from LTBI ("Disease"),
#' TB from some combination of other diseases and possibly LTBI
#' ("Disease/Other Diseases"), TB from Human Immunodeficiency Virus
#' ("Disease/HIV"), TB from pneumonia ("Disease/Pneumonia"),
#' or identify risk of progression to TB ("risk"), risk of TB treatment
#' failure ("failure"), or classify treatment responses
#' (i.e., failures from cures, "response"). Required.
#' @param tissuetype a \code{character} string that denotes whether the
#' signature was developed using samples of either whole blood/paxgene
#' ("whole blood") or peripheral blood mononuclear cells ("PBMC").
#' Due to the manipulation of cells inherently required to obtain PBMCs,
#' many scientists prefer to use only whole blood samples for analysis.
#' Accepts "whole blood", "PBMC" or "mixed". Required.
#' @param saveobjs \code{logical}. If \code{TRUE}, the contents of the data
#' file (\code{TBsignatures}, \code{TBcommon}, \code{sigAnnotData},
#' \code{common_sigAnnotData}) will be overwritten and updated to include
#' the new signature. If \code{FALSE}, no files will be overwritten, but
#' you can check function output for errors before writing RDS objects
#' by setting \code{views = TRUE}. Default is \code{saveobjs = FALSE}.
#' @param views logical. If \code{TRUE}, all objects will be sent to a
#' data view in a new window to check for errors. Default is \code{TRUE}.
#'
#' @return Either data objects \code{TBsignatures}, \code{TBcommon},
#' \code{sigAnnotData}, and \code{common_sigAnnotData} will be updated with the
#' new signature and overwritten if \code{saveobjs = FALSE}, or no output will
#' be produced except errors and messages for checking that the function
#' runs correctly given the inputs.
#'
#' @examples
#' # Mock example signature
#' TBSignatureProfiler:::addTBsignature(sigsymbols = c("GBP5", "BATF2", "GZMA"),
#'                             authname = "Odom",
#'                             signame_common = NULL,
#'                             sigtype = "Disease/HIV",
#'                             tissuetype = "PBMC",
#'                             saveobjs = FALSE,
#'                             views = FALSE)

addTBsignature <- function(sigsymbols, authname,
                         signame_common = NULL, sigtype, tissuetype,
                         saveobjs = FALSE, views = TRUE) {
  # Checks and renaming
  if (!is.character(sigsymbols)) {
    stop("'sigsymbols' must be a character vector.")
  } else if (!is.character(authname)) {
    stop("'authname' must be a character string.")
  } else if (length(sigsymbols) <= 1) {
    stop("Currently the TBSP only accepts gene signatures with 2 or more ",
               "gene symbols.")
  } else if (!is.logical(saveobjs)) {
    stop("'saveobjs' must be logical.")
  }

  diseasetypes <- c("Disease", "Disease/HIV", "Disease/Other Diseases",
                    "Disease/Pneumonia", "failure", "response", "risk")
  replacetypes <- c("Disease", "HIV", "OD", "PNA", "FAIL", "RES", "RISK")
  if (!(sigtype %in% diseasetypes)) {
    stop("Input string for 'sigtype' not recognized. You must denote ",
    "'sigtype' as one of the following: 'Disease', 'Disease/HIV', ",
    "'Disease/Other Diseases', 'Disease/Pneumonia', 'failure', ",
    "'response', 'risk'")
  } else {
    sigtype <- replacetypes[which(sigtype == diseasetypes)]
  }

  alltissues <- c("mixed", "PBMC", "whole blood")
  if (!(tissuetype %in% alltissues)) {
    stop("'tissuetype' must be one of the following: ",
         "'mixed', 'PBMC', 'whole blood'")
  } else if (tissuetype == "whole blood") tissuetype <- "whole"

  # Create signature name
  ngenes <- length(sigsymbols)
  if (sigtype == "Disease") {
    signame_TBSP <- paste(authname, ngenes, sep = "_")
  } else signame_TBSP <- paste(authname, sigtype, ngenes, sep = "_")

  if (nchar(signame_TBSP) > 17) {
    howmanyletters <- 17 - (2 + nchar(as.character(ngenes)) + length(sigtype))
    stop("Resulting signature name is too long. ",
         "Please shorten or abbreviate authname to ", howmanyletters,
         " characters.")
  } else if (signame_TBSP %in% names(TBsignatures)) {
    stop("There is already a signature with that name. ",
         "Please alter the authname to distinguish it from the previous ",
         "signature.")
  } else message("The assigned signature name is ", signame_TBSP)

  if (is.null(signame_common)) {
    message("No alternative signature name was provided")
    signame_common <- signame_TBSP
  } else message("Note: 'signame_common' differs from 'signame_TBSP'")

  # Add signature to TBsignatures
  TBsignatures$`NEWSIG` <-  sigsymbols
  names(TBsignatures)[names(TBsignatures) == "NEWSIG"] <- signame_TBSP
  # Alphabetize
  index <- order(names(TBsignatures))
  new_names <- names(TBsignatures)[index]
  TBsignatures <- TBsignatures[new_names]
  # Save new object
  if (saveobjs) save(TBsignatures, file = "data/TBsignatures.rda")
  if (views) utils::View(names(TBsignatures))
  message("TBsignatures object updated")

  # Add signature to TBcommon
  TBcommon$`NEWSIG` <-  sigsymbols
  names(TBcommon)[names(TBcommon) == "NEWSIG"] <- signame_common
  # reorder names
  new_names <- names(TBcommon)[index]
  TBcommon <- TBcommon[new_names]
  if (saveobjs) save(TBcommon, file = "data/TBcommon.rda")
  if (views) utils::View(names(TBcommon))
  message("TBcommon object updated")

  # Add to annotation to sigAnnotData
  newrow <- c("names" = signame_TBSP, "disease" = sigtype,
              "tissue_type" = tissuetype)
  sigAnnotData <- as.data.frame(rbind(sigAnnotData, newrow))[index, ]
  sigAnnotData$disease <- as.factor(sigAnnotData$disease)
  sigAnnotData$tissue_type <- as.factor(sigAnnotData$tissue_type)
  rownames(sigAnnotData) <- seq_len(nrow(sigAnnotData))
  if (all(sigAnnotData$names != names(TBsignatures))) {
    stop("The names in the sigAnnotData do not match the order ",
         "and/or names of the TBsignatures object.")
  } else if (saveobjs) save(sigAnnotData, file = "data/sigAnnotData.rda")
  if (views) utils::View(sigAnnotData)
  message("sigAnnotData object updated")

  # Add annotation to sigAnnotData_common
  newrow <- c("names" = signame_common, "disease" = sigtype,
              "tissue_type" = tissuetype)
  common_sigAnnotData <- as.data.frame(rbind(common_sigAnnotData,
                                             newrow))[index, ]
  common_sigAnnotData$disease <- as.factor(common_sigAnnotData$disease)
  common_sigAnnotData$tissue_type <- as.factor(
    common_sigAnnotData$tissue_type)
  rownames(common_sigAnnotData) <- seq_len(nrow(common_sigAnnotData))
  if (all(common_sigAnnotData$names != names(TBcommon))) {
    stop("The names in the common_sigAnnotData table do not match ",
         "the order and/or names of the TBcommon object.")
  } else if (saveobjs) save(common_sigAnnotData,
                            file = "data/common_sigAnnotData.rda")
  if (views) utils::View(common_sigAnnotData)
  message("common_sigAnnotData updated")

  if (saveobjs) message("All objects updated and saved to data folder.")
}


