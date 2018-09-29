#TODO: input1: signature data, annotation data frame (matrix of signature data and matrix of annotations)
#TODO: input2: SummarizedExperiment object, list of columns of signature data, list of columns of annotation data
#TODO: Better example
#TODO: try with one signature
# signatureHeatmap(india_relapse_GSVA, signatureColNames = colnames(colData(india_relapse_GSVA))[71:85],
# annotationColNames=c("visit", "subjtype"))
#' Plot a Heatmap of Signature Genes
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
#' mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
#'                              dimnames = list(TBsignatures$Predict29,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10))))
#' res <- runTBsigProfiler(mat_testdata, algorithm = "GSVA", parallel.sz = 1)
signatureHeatmap <- function(inputData, annotationData, name="Signatures",
                             signatureColNames, annotationColNames, colList,
                             scale=FALSE, showColumnNames=TRUE, showRowNames=TRUE,
                             colorSets=c("Set1", "Set2", "Set3", "Pastel1",
                                         "Pastel2", "Accent", "Dark2",
                                         "Paired")) {
  colList <- list()
  colorSetNum <- 1
  for (annot in annotationColNames){
    condLevels <- unique(SummarizedExperiment::colData(inputData)[, annot])
    if (length(condLevels) > 8){
      colors <- distinctColors(length(condLevels))
    } else {
      colors <- RColorBrewer::brewer.pal(8, colorSets[colorSetNum])
      colorSetNum <- colorSetNum + 1
    }
    colList[[annot]] <- stats::setNames(colors[seq_along(condLevels)],
                                        condLevels)
  }

  topha2 <- ComplexHeatmap::HeatmapAnnotation(
    df = data.frame(SummarizedExperiment::colData(inputData)[, annotationColNames]),
    col = colList, height = grid::unit(0.4 * length(annotationColNames), "cm"),
    show_legend = TRUE, show_annotation_name = TRUE)
  sigresults <- t(as.matrix(SummarizedExperiment::colData(inputData)[, signatureColNames]))
  keyname <- "Score"
  if (scale){
    sigresults <- t(scale(t(sigresults)))
    keyname <- "Scaled\nScore"
  }
  ComplexHeatmap::draw(
    ComplexHeatmap::Heatmap(sigresults, column_title = name,
                            show_column_names = showColumnNames,
                            show_row_names = showRowNames,
                            top_annotation = topha2, name = keyname),
    annotation_legend_side = "bottom"
  )
}

#TODO: Better example
#TODO: test with 1 input
#TODO: support other input method
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
#'
#' @return A ggplot2 boxplot of the signature data using the annotation
#' provided.
#'
#' @export
#'
#' @examples
#' mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
#'                              dimnames = list(TBsignatures$Predict29,
#'                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                              paste0("sample", 1:10))))
#' res <- runTBsigProfiler(mat_testdata, algorithm = "GSVA", parallel.sz = 1)
signatureBoxplot <- function(inputData, annotationData, signatureColNames,
                             annotationColName, name="Signatures", scale=FALSE) {
  pathwaydata <- t(as.matrix(SummarizedExperiment::colData(inputData)[, signatureColNames]))
  if (scale) {
    pathwaydata <- t(scale(t(pathwaydata)))
  }
  boxplotdf <- data.frame(t(pathwaydata),
                          Group = SummarizedExperiment::colData(inputData)[, annotationColName])
  boxplotdfm <- reshape2::melt(boxplotdf, value.name = "Score",
                               variable.name = "Signature", id.vars = "Group")
  ggplot2::ggplot(boxplotdfm, ggplot2::aes_string("Group", "Score")) +
    ggplot2::facet_wrap(~Signature, scales = 'free') +
    ggplot2::geom_boxplot(outlier.shape = NA, ggplot2::aes_string(fill = "Group")) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1)) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::ggtitle(name)
}

signatureGeneHeatmap <- function(inputData, geneList, annotationData,
                                 signatureColNames, annotationColNames) {

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
