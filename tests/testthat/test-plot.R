context("test-plot")

set.seed(1234)
annotdata <- data.frame(sample = c(rep("down", 5), rep("up", 5)),
                        samplename = paste0("sample", 1:10),
                        row.names = paste0("sample", 1:10))
mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
                             dimnames = list(TBsignatures$Zak_RISK_16,
                                             paste0("sample", 1:10))),
                      matrix(rnorm(1000), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
df_testdata <- data.frame(mat_testdata)
SE_testdata <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = mat_testdata),
  colData = S4Vectors::DataFrame(sample = c(rep("down", 5), rep("up", 5)),
                                 samplename = paste0("sample", 1:10),
                                 samplenum = c(1:5, 1:5)))
SE_wres <- runTBsigProfiler(SE_testdata, useAssay = "data",
                            signatures = list(sig1 = paste0("gene", 1:10)),
                            algorithm = "ssGSEA", parallel.sz = 1)
dfres <- runTBsigProfiler(SE_testdata, useAssay = "data",
                          outputFormat = "data.frame",
                          signatures = list(sig1 = paste0("gene", 1:10)),
                          algorithm = "ssGSEA", parallel.sz = 1)
df_errornames <- dfres
colnames(df_errornames)[1] <- "errorname"
terrorname <- t(df_errornames)

test_that("Missing annotationData", {
  #signatureHeatmap
  expect_error(
    signatureHeatmap(inputData = SE_wres),
    "argument \"signatureColNames\" is missing, with no default")
  expect_error(
    signatureHeatmap(inputData = SE_wres, signatureColNames = "sig1"),
    "argument \"annotationColNames\" is missing, with no default")
  expect_error(
    signatureHeatmap(inputData = SE_wres, signatureColNames = "sig1",
                     annotationColNames = "group"),
    "Annotation column name not found in inputData.")
  expect_error(
    signatureHeatmap(inputData = SE_wres, signatureColNames = c("sig1", "sig1"),
                     annotationColNames = "group"),
    "Duplicate signature column name is not supported.")
  expect_error(
    signatureHeatmap(inputData = SE_wres, signatureColNames = c("sig2"),
                     annotationColNames = "group"),
    "Signature column name not found in inputData.")
  expect_error(
    signatureHeatmap(inputData = df_errornames, annotationData = annotdata),
    "Annotation data and signature data does not match.")
  expect_error(
    signatureHeatmap(inputData = terrorname, annotationData = annotdata),
    "Annotation data and signature data does not match.")
  expect_error(
    signatureHeatmap(inputData = df_errornames[, 1:5],
                     annotationData = annotdata),
    "Annotation data and signature data does not match.")

  #signatureBoxplot
  expect_error(
    signatureBoxplot(inputData = SE_wres),
    "argument \"signatureColNames\" is missing, with no default")
  expect_error(
    signatureBoxplot(inputData = SE_wres,
                     signatureColNames = "sig1"),
    "argument \"annotationColName\" is missing, with no default")
  expect_error(
    signatureBoxplot(inputData = SE_wres,
                     annotationColName = "sample",
                     signatureColNames = "sig2"),
    "Signature column name not found in inputData.")
  expect_error(
    signatureBoxplot(inputData = SE_wres,
                     annotationColName = "sample",
                     signatureColNames = c("sig1", "sig1")),
    "Duplicate signature column name is not supported.")
  expect_error(
    signatureBoxplot(inputData = SE_wres,
                     annotationColName = c("sample", "sample"),
                     signatureColNames = "sig1"),
    "You must specify a single annotation column name to color boxplots by.")
  expect_error(
    signatureBoxplot(inputData = SE_wres,
                     annotationColName = "samplea",
                     signatureColNames = "sig1"),
    "Annotation column name not found in inputData.")
  expect_error(
    signatureBoxplot(inputData = dfres,
                     annotationData = annotdata),
    "annotationData must have only one column.")
  expect_error(
    signatureBoxplot(inputData = df_errornames,
                     annotationData = annotdata[, 1, drop = FALSE]),
    "Annotation data and signature data does not match.")
  expect_error(
    signatureBoxplot(inputData = terrorname,
                     annotationData = annotdata[, 1, drop = FALSE]),
    "Annotation data and signature data does not match.")
  expect_error(
    signatureBoxplot(inputData = SE_wres, annotationColName = "samplename",
                     signatureColNames = "sig1"),
    "Too many levels in the annotation data. The boxplot can contain a maximum of 9 levels")
  expect_error(
    signatureBoxplot(inputData = df_errornames[, 1:5],
                     annotationData = annotdata[, 1, drop = FALSE]),
    "Annotation data and signature data does not match.")
  expect_error(
    distinctColors(10, hues = "peuce"),
    "Only color names listed in the \'color\' function can be used in \'hues\'")

  #signatureGeneHeatmap
  expect_error(
    signatureGeneHeatmap(),
    "argument \"inputData\" is missing, with no default")
  expect_error(
    signatureGeneHeatmap(inputData = SE_wres),
    "argument \"useAssay\" is missing, with no default")
  expect_error(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data"),
    "argument \"sigGenes\" is missing, with no default")
  expect_error(
    signatureGeneHeatmap(inputData = SE_wres, signatureColNames = "sig1",
                         annotationColNames = "sample"),
    "argument \"useAssay\" is missing, with no default")
  expect_error(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         signatureColNames = "sig1",
                         annotationColNames = "sample"),
    "argument \"sigGenes\" is missing, with no default")
  expect_error(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         sigGenes = TBsignatures$Zak_RISK_16,
                         annotationColNames = c("sample"),
                         colList = list("wrong" = c("yes" = "blue"))),
    "The colList is out of sync with the annotation columns")
})

test_that("SummarizedExperiment Plot Works", {
  expect_is(
    signatureHeatmap(inputData = SE_wres, signatureColNames = "sig1",
                     annotationColNames = "sample"),
    "HeatmapList"
  )
  expect_is(
    signatureHeatmap(inputData = SE_wres, signatureColNames = "sig1",
                     annotationColNames = "sample", scale = TRUE),
    "HeatmapList"
  )

  #signatureBoxplot
  expect_is(
    signatureBoxplot(inputData = SE_wres, annotationColName = "sample",
                     signatureColNames = "sig1"),
    "ggplot"
  )
  expect_is(
    signatureBoxplot(inputData = SE_wres, annotationColName = "sample",
                     signatureColNames = "sig1", scale = TRUE),
    "ggplot"
  )
  expect_is(
    signatureBoxplot(inputData = SE_wres, annotationColName = "samplenum",
                     signatureColNames = "sig1"),
    "ggplot"
  )

  expect_is(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         sigGenes = TBsignatures$Zak_RISK_16),
    "HeatmapList"
  )
  expect_is(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         sigGenes = TBsignatures$Zak_RISK_16,
                         signatureColNames = "sig1"),
    "HeatmapList"
  )
  expect_is(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         sigGenes = TBsignatures$Zak_RISK_16,
                         annotationColNames = c("sample")),
    "HeatmapList"
  )
  expect_is(
    signatureGeneHeatmap(inputData = SE_wres, useAssay = "data",
                         sigGenes = TBsignatures$Zak_RISK_16,
                         signatureColNames = "sig1",
                         annotationColNames = c("sample", "samplename")),
    "HeatmapList"
  )
})

test_that("DataFrame Plot Works", {
  expect_is(
    signatureHeatmap(inputData = dfres,
                     annotationData = annotdata),
    "HeatmapList"
  )
  expect_is(
    signatureBoxplot(inputData = dfres,
                     annotationData = annotdata[, 1, drop = FALSE]),
    "ggplot"
  )
})
