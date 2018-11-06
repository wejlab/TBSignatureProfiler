#' Plot a Heatmap of Signature Scores
#'
#' @param inputData Either a SummarizedExperiment object that contains the
#' signature data and annotation as colData columns, or a data.frame or matrix
#' of signature data. Required.
#' @param annotationData If inputData is a data.frame or matrix of signature
#' data, a data.frame or matrix of annotation data.
#' @param name The name of the heatmap. The default is "Signatures".
#' @param signatureColNames If inputData is a SummarizedExperiment, the column
#' names in colData that contain the signature data.
#' @param annotationColNames If inputData is a SummarizedExperiment, the column
#' names in colData that contain the annotation data.
#' @param colList Custom color information as a named list to pass to
#' ComplexHeatmap. By default, ColorBrewer color sets will be used. See
#' colorSets for additional options.
#' @param scale Scale the signature data. The default is FALSE.
#' @param showColumnNames Show columns names (sample names) in the heatmap. The
#' default is TRUE.
#' @param showRowNames Show row names(signature names) in the heatmap. The
#' default is TRUE.
#' @param colorSets By default, this function will use the color sets in the
#' order listed below for annotation information. Replace this with the sets
#' in order that you want to use them, or provide custom color sets with
#' colList.
#'
#' @return A ComplexHeatmap plot
#'
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' #generate some artificial data that shows a difference in ACS_COR_16
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$ACS_COR_16,
#'                                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                                              paste0("sample", 1:10))))
#' #Create a SummarizedExperiment object that contains the data
#' testdataSE <- SummarizedExperiment(assays=SimpleList(data = mat_testdata),
#'                                    colData=DataFrame(sample = c(rep("down", 5),
#'                                                                 rep("up", 5))))
#' #Run profiler using GSVA and ssGSEA on ACS_COR_16
#' res <- runTBsigProfiler(testdataSE, useAssay = "data",
#'                         signatures = TBsignatures["ACS_COR_16"],
#'                         algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1,
#'                         combineSigAndAlgorithm = TRUE)
#' #Plot heatmap of the results
#' signatureHeatmap(res, signatureColNames = c("GSVA_ACS_COR_16",
#'                                             "ssGSEA_ACS_COR_16"),
#'                  annotationColNames = "sample", scale=TRUE,
#'                  showColumnNames = FALSE)
signatureHeatmap <- function(inputData, annotationData, name="Signatures",
                             signatureColNames, annotationColNames,
                             colList=list(), scale=FALSE, showColumnNames=TRUE,
                             showRowNames=TRUE, colorSets=c("Set1", "Set2",
                             "Set3", "Pastel1", "Pastel2", "Accent", "Dark2",
                             "Paired")) {
  if (methods::is(inputData, "SummarizedExperiment")){
    if (any(duplicated(signatureColNames))){
      stop("Duplicate signature column name is not supported.")
    }
    if (!all(signatureColNames %in% colnames(SummarizedExperiment::colData(inputData)))){
      stop("Signature column name not found in inputData.")
    }
    if (!all(annotationColNames %in% colnames(SummarizedExperiment::colData(inputData)))){
      stop("Annotation column name not found in inputData.")
    }
    annotationData <- SummarizedExperiment::colData(inputData)[, annotationColNames, drop = FALSE]
    inputData <-  SummarizedExperiment::colData(inputData)[, signatureColNames, drop = FALSE]
  } else {
    annotationColNames <- colnames(annotationData)
  }
  #if number of rows equal number of rows
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

  if (length(colList) == 0){
    colorSetNum <- 1
    for (annot in annotationColNames){
      if (is.numeric(annotationData[, annot])){
        t1min <- min(annotationData[, annot], na.rm = TRUE)
        t1max <- max(annotationData[, annot], na.rm = TRUE)
        colList[[annot]] <- circlize::colorRamp2(c(t1min, t1max),
                                                 c("white", "blue"))
      } else {
        condLevels <- unique(annotationData[, annot][!is.na(annotationData[, annot])])
        if (length(condLevels) > 8){
          colors <- distinctColors(length(condLevels))
        } else {
          colors <- RColorBrewer::brewer.pal(8, colorSets[colorSetNum])
          colorSetNum <- colorSetNum + 1
        }
        colList[[annot]] <- stats::setNames(colors[seq_along(condLevels)],
                                            condLevels)
      }
    }
  }

  topha2 <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(annotationData),
    col = colList, height = grid::unit(0.4 * length(annotationColNames), "cm"),
    show_legend = TRUE, show_annotation_name = TRUE)
  sigresults <- t(as.matrix(inputData))
  keyname <- "Score"
  if (scale){
    sigresults <- t(scale(t(sigresults)))
    keyname <- "Scaled\nScore"
  }
  return(ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(sigresults, column_title = name,
                            show_column_names = showColumnNames,
                            show_row_names = showRowNames,
                            top_annotation = topha2, name = keyname),
    annotation_legend_side = "bottom"
  ))
}

#' Plot a Boxplot of Signature Genes
#'
#' @param inputData Either a SummarizedExperiment object that contains the
#' signature data and annotation as colData columns, or a data.frame or matrix
#' of signature data. Required.
#' @param annotationData If inputData is a data.frame or matrix of signature
#' data, a data.frame or matrix of annotation data.
#' @param signatureColNames If inputData is a SummarizedExperiment, the column
#' names in colData that contain the signature data.
#' @param annotationColName If inputData is a SummarizedExperiment, the column
#' name in colData to use to make the boxplot.
#' @param name The name of the boxplot. The default is "Signatures".
#' @param scale Scale the signature data. The default is FALSE.
#' @param includePoints Include points over the boxplots. The default is TRUE.
#'
#' @return A ggplot2 boxplot of the signature data using the annotation
#' provided.
#'
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' #generate some artificial data that shows a difference in ACS_COR_16
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$ACS_COR_16,
#'                                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                                              paste0("sample", 1:10))))
#' #Create a SummarizedExperiment object that contains the data
#' testdataSE <- SummarizedExperiment(assays=SimpleList(data = mat_testdata),
#'                                    colData=DataFrame(sample = c(rep("down", 5),
#'                                                                 rep("up", 5))))
#' #Run profiler using GSVA and ssGSEA on ACS_COR_16
#' res <- runTBsigProfiler(testdataSE, useAssay = "data",
#'                         signatures = TBsignatures["ACS_COR_16"],
#'                         algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1,
#'                         combineSigAndAlgorithm = TRUE)
#' #Plot boxplot of the results
#' signatureBoxplot(res, signatureColNames = c("GSVA_ACS_COR_16",
#'                                             "ssGSEA_ACS_COR_16"),
#'                  annotationColName = "sample", name = "ACS_COR_16 Signatures")
signatureBoxplot <- function(inputData, annotationData, signatureColNames,
                             annotationColName, name="Signatures", scale=FALSE,
                             includePoints=TRUE) {
  if (methods::is(inputData, "SummarizedExperiment")){
    if (any(duplicated(signatureColNames))){
      stop("Duplicate signature column name is not supported.")
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
  if (length(levels(annotationData[, 1])) > 9){
    stop("Too many levels in the annotation data. The boxplot can contain a maximum of 9 levels")
  }
  #if number of rows equal number of rows
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
  pathwaydata <- t(inputData)
  if (scale) {
    pathwaydata <- t(scale(t(pathwaydata)))
  }
  boxplotdf <- data.frame(t(pathwaydata),
                          Group = annotationData[, 1])
  boxplotdfm <- reshape2::melt(boxplotdf, value.name = "Score",
                               variable.name = "Signature", id.vars = "Group")
  theplot <- ggplot2::ggplot(boxplotdfm, ggplot2::aes_string("Group", "Score")) +
    ggplot2::facet_wrap(~Signature, scales = 'free') +
    ggplot2::geom_boxplot(outlier.shape = NA, ggplot2::aes_string(fill = "Group")) +
    ggplot2::theme_classic()
  if (includePoints) {
    theplot <- theplot + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1))
  }
  return(theplot +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::ggtitle(name))
}

#' Plot a Heatmap of a single signature scores and gene expression
#'
#' @param inputData A SummarizedExperiment object that contains the
#' signature data and annotation as colData columns. Required.
#' @param useAssay The assay to use for the gene expression data. Required.
#' @param sigGenes The genes in the signature to use in the heatmap. For inbuilt
#' signatures, you can use TBsignatures, (ex: TBsignatures[["ACS_COR"]]).
#' Required.
#' @param name The name of the heatmap. The default is "Signature".
#' @param signatureColNames The column names in colData that contain the
#' signature data.
#' @param annotationColNames The column names in colData that contain the
#' annotation data.
#' @param scale Scale the gene expression data. The default is TRUE.
#' @param showColumnNames Show columns names (sample names) in the heatmap. The
#' default is TRUE.
#' @param showRowNames Show row names(signature names) in the heatmap. The
#' default is TRUE.
#' @param colList Custom color information for annotation data as a named list
#' to pass to ComplexHeatmap. By default, ColorBrewer color sets will be used.
#' See colorSets for additional options.
#' @param colorSets By default, this function will use the color sets in the
#' order listed below for annotation information. Replace this with the sets
#' in order that you want to use them, or provide custom color sets with
#' colList.
#'
#' @return  A ComplexHeatmap plot
#'
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' #generate some artificial data that shows a difference in ACS_COR_16
#' mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
#'                              dimnames = list(TBsignatures$ACS_COR_16,
#'                                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                                              paste0("sample", 1:10))))
#' #Create a SummarizedExperiment object that contains the data
#' testdataSE <- SummarizedExperiment(assays=SimpleList(data = mat_testdata),
#'                                    colData=DataFrame(sample = c(rep("down", 5),
#'                                                                 rep("up", 5))))
#' #Run profiler using GSVA and ssGSEA on ACS_COR_16
#' res <- runTBsigProfiler(testdataSE, useAssay = "data",
#'                         signatures = TBsignatures["ACS_COR_16"],
#'                         algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1,
#'                         combineSigAndAlgorithm = TRUE)
#' #Plot a heatmap of signature genes and pathway predictions
#' signatureGeneHeatmap(res, useAssay = "data",
#'                      sigGenes = TBsignatures[["ACS_COR_16"]],
#'                      signatureColNames = c("GSVA_ACS_COR_16",
#'                                            "ssGSEA_ACS_COR_16"),
#'                      annotationColNames = c("sample"), showColumnNames = FALSE,
#'                      name = "ACS_COR_16")
signatureGeneHeatmap <- function(inputData, useAssay, sigGenes,
                                 name="Signature", signatureColNames=NULL,
                                 annotationColNames = NULL, scale=TRUE,
                                 showColumnNames=TRUE, showRowNames=TRUE,
                                 colList=list(), colorSets=c("Set1", "Set2",
                                 "Set3", "Pastel1", "Pastel2", "Accent",
                                 "Dark2", "Paired")) {
  if (!is.null(signatureColNames)){
    pathwaycols <- list()
    pathwaydata <- data.frame(SummarizedExperiment::colData(inputData)[, signatureColNames, drop = FALSE])
    for (i in colnames(pathwaydata)){
      t1min <- min(pathwaydata[, i], na.rm = TRUE)
      t1max <- max(pathwaydata[, i], na.rm = TRUE)
      pathwaycols[[i]] <- circlize::colorRamp2(c(t1min, mean(c(t1min, t1max)), t1max),
                                               c("darkgreen", "white", "darkorange"))
    }
  } else {
    pathwaycols <- NULL
    pathwaydata <- NULL
  }

  annotationData <- data.frame(SummarizedExperiment::colData(inputData)[, annotationColNames, drop = FALSE])

  if (!is.null(annotationColNames)) {
    if (length(colList) == 0){
      colorSetNum <- 1
      for (annot in annotationColNames){
        if (is.numeric(annotationData[, annot])){
          t1min <- min(annotationData[, annot], na.rm = TRUE)
          t1max <- max(annotationData[, annot], na.rm = TRUE)
          colList[[annot]] <- circlize::colorRamp2(c(t1min, t1max),
                                                   c("white", "blue"))
        } else {
          condLevels <- unique(annotationData[, annot][!is.na(annotationData[, annot])])
          if (length(condLevels) > 8){
            colors <- distinctColors(length(condLevels))
          } else {
            colors <- RColorBrewer::brewer.pal(8, colorSets[colorSetNum])
            colorSetNum <- colorSetNum + 1
          }
          colList[[annot]] <- stats::setNames(colors[seq_along(condLevels)],
                                              condLevels)
        }
      }
    }
    if (any(annotationColNames != names(colList))) {
      stop("The colList is out of sync with the annotation columns")
    }
  } else {
    colList <- NULL
  }

  if (!is.null(pathwaydata) | !is.null(annotationColNames)) {
    if (!is.null(annotationColNames) & !is.null(pathwaydata)){
      annotDF <- cbind(data.frame(SummarizedExperiment::colData(inputData)[, annotationColNames, drop = FALSE]),
                       pathwaydata)
    } else if (!is.null(annotationColNames)) {
      annotDF <- cbind(data.frame(SummarizedExperiment::colData(inputData)[, annotationColNames, drop = FALSE]))
    } else {
      annotDF <- pathwaydata
    }
    topha <- ComplexHeatmap::HeatmapAnnotation(
      df = annotDF,
      col = c(colList, pathwaycols),
      height = grid::unit(0.4 * length(c(annotationColNames, signatureColNames)), "cm"), show_legend = TRUE, show_annotation_name = TRUE)
  } else {
    topha <- NULL
  }
  heatdata <- SummarizedExperiment::assay(inputData, useAssay)[sigGenes[sigGenes %in% rownames(inputData)], ]
  heatname <- useAssay
  if (scale){
    heatdata <- heatdata[rowSums(heatdata) != 0, ]
    heatdata <- t(scale(t(heatdata)))
    heatname <- paste("Scaled", heatname, sep = "\n")
  }
  return(ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(
      heatdata, show_column_names = showColumnNames,
      show_row_names = showRowNames, top_annotation = topha,
      name = heatname, column_title = name),
    annotation_legend_side = "bottom")
  )
}

#' Generate a distinct palette for coloring different clusters
#'
#' @param n Integer; Number of colors to generate
#' @param hues Character vector of R colors available from the colors()
#' function. These will be used as the base colors for the clustering scheme.
#' Different saturations and values (i.e. darkness) will be generated for each
#' hue.
#' @param saturation.range Numeric vector of length 2 with values between 0 and
#' 1. Default: c(0.25, 1)
#' @param value.range Numeric vector of length 2 with values between 0 and 1.
#' Default: c(0.5, 1)
#' @return A vector of distinct colors that have been converted to  HEX from
#' HSV.
#' @export
distinctColors <- function(n, hues = c("red", "cyan", "orange", "blue",
                                       "yellow", "purple", "green", "magenta"),
                           saturation.range = c(0.7, 1),
                           value.range = c(0.7, 1)) {
  #Adapted from compbiomed/celda, thanks to all celda developers
  if (!(all(hues %in% grDevices::colors()))) {
    stop("Only color names listed in the 'color'",
         " function can be used in 'hues'")
  }

  ## Convert R colors to RGB and then to HSV color format
  hues.hsv <- grDevices::rgb2hsv(grDevices::col2rgb(hues))

  ## Calculate all combination of saturation/value pairs
  ## Note that low saturation with low value (i.e. high darkness) is too dark
  ## for all hues
  ## Likewise, high saturation with high value (i.e. low darkness) is hard to
  ## distinguish
  ## Therefore, saturation and value are set to be anticorrelated
  num.vs <- ceiling(n / length(hues))
  s <- seq(from = saturation.range[1], to = saturation.range[2],
           length = num.vs)
  v <- seq(from = value.range[2], to = value.range[1], length = num.vs)

  ## Create all combination of hues with saturation/value pairs
  new.hsv <- c()
  for (i in 1:num.vs) {
    temp <- rbind(hues.hsv[1, ], s[i], v[i])
    new.hsv <- cbind(new.hsv, temp)
  }

  ## Convert to hex
  col <- grDevices::hsv(new.hsv[1, ], new.hsv[2, ], new.hsv[3, ])

  return(col[1:n])
}
