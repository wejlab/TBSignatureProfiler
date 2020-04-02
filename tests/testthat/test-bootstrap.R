context("test-bootstrap")

set.seed(1234)
mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
                             dimnames = list(TBsignatures$Zak_RISK_16,
                                             paste0("sample", 1:10))),
                      matrix(rnorm(1000), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
df_testdata <- data.frame(mat_testdata)
SEtestdata <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = mat_testdata))
SummarizedExperiment::colData(SEtestdata)[, 1] <- factor(c(rep("yes", 5),
                                                           rep("no", 5)))
SummarizedExperiment::colData(SEtestdata)[, 2] <- c(rep("yes", 5),
                                                    rep("no", 5))
names(SummarizedExperiment::colData(SEtestdata)) <- c("Disease", "Notfac")
scored_testdata <- runTBsigProfiler(SEtestdata, useAssay = "data",
                                    signatures = list(sig1 = paste0("gene", 1:10)),
                                    algorithm = "ssGSEA", parallel.sz = 1)
annot_data <- data.frame(factor(c(rep("yes", 5), rep("no", 5))))
colnames(annot_data) <- "Disease"
rownames(annot_data) <- paste0("sample", 1:10)
scored_df <- runTBsigProfiler(df_testdata, output = "data.frame",
                              signatures = list(sig1 = paste0("gene", 1:10)),
                              algorithm = "ssGSEA", parallel.sz = 1)


test_that("Run bootstrapAUC", {
  expect_output(
    bootstrapAUC(SE_scored = scored_testdata,
                 signatureColNames = "sig1",
                 annotationColName = "Disease",
                 num.boot = 5)
  )
  expect_output(
    bootstrapAUC(SE_scored = scored_testdata,
                 signatureColNames = "sig1",
                 annotationColName = "Notfac",
                 num.boot = 5)
  )
})

test_that("Run tableAUC", {
  expect_output(
    tableAUC(SE_scored = scored_testdata,
                 signatureColNames = "sig1",
                 annotationColName = "Disease",
                 num.boot = 5)
  )
  expect_is(
    tableAUC(SE_scored = scored_testdata,
                 signatureColNames = "sig1",
                 annotationColName = "Disease",
                 num.boot = 5,
                 output = "DataTable"),
    "datatables"
  )
  expect_is(
    tableAUC(SE_scored = scored_testdata,
             signatureColNames = "sig1",
             annotationColName = "Disease",
             num.boot = 5,
             output = "data.frame"),
    "data.frame"
  )
})

test_that("compareBoxplots", {
  expect_output(
    compareBoxplots(SE_scored = scored_testdata,
                                signatureColNames = "sig1",
                                annotationColName = "Disease",
                                num.boot = 5)
    )
  expect_output(
    compareBoxplots(SE_scored = scored_testdata,
                    signatureColNames = "sig1",
                    annotationColName = "Disease",
                    num.boot = 5,
                    rotateLabels = TRUE)
  )
})

test_that("Run signatureROCplot", {
  expect_is(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Disease"),
    "gg"
  )
  expect_is(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Disease",
                     scale = TRUE),
    "gg"
  )
  expect_is(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = c("sig1", "sig1"),
                     annotationColName = "Disease"),
    "gg"
  )
  expect_error(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = "sig2",
                     annotationColName = "Disease"),
    "Signature column name not found in inputData."
  )
  expect_error(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Not_here"),
    "Annotation column name not found in inputData."
  )
  expect_error(
    signatureROCplot(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = c("Disease", "Notfac")),
    "You must specify a single annotation column name with which
    to color boxplots."
  )
})

test_that("Run signatureROCplot_CI", {
  expect_is(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Disease"),
    "gg"
  )
  expect_is(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Disease",
                     scale = TRUE),
    "gg"
  )
  expect_is(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = c("sig1", "sig1"),
                     annotationColName = "Disease"),
    "gg"
  )
  expect_error(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = "sig2",
                     annotationColName = "Disease"),
    "Signature column name not found in inputData."
  )
  expect_error(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = "Not_here"),
    "Annotation column name not found in inputData."
  )
  expect_error(
    signatureROCplot_CI(inputData = scored_testdata,
                     signatureColNames = "sig1",
                     annotationColName = c("Disease", "Notfac")),
    "You must specify a single annotation column name with which
    to color boxplots."
  )
})

Mat_testdata <- rbind(matrix(c(rpois(80, lambda = 50), 
                               rpois(80, lambda = 50) + 5), 16, 10,
                             dimnames = list(TBsignatures$Zak_RISK_16,
                                             paste0("sample", 1:10))),
                      matrix(rpois(1000, lambda = 50), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
SEtestdata2 <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = Mat_testdata))

test_that("Run mk_Assay", {
  expect_is(
    mkAssay(SEtestdata2, input_name = "data", log = TRUE,
            counts_to_CPM = TRUE),
    "SummarizedExperiment"
  )
  expect_error(
    mkAssay(SEtestdata2, input_name = "data", log = FALSE,
            counts_to_CPM = FALSE),
    "At least counts_to_CPM or log must be TRUE."
  )
  expect_error(
    mkAssay(SEtestdata2, input_name = "not_here"),
    "input_name must be an SE_obj assay"
  )
  
})