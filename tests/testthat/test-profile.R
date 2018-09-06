context("test-profile")

set.seed(1234)
mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
                             dimnames = list(TBsignatures$Predict29,
                                             paste0("sample", 1:10))),
                      matrix(rnorm(1000), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
df_testdata <- data.frame(mat_testdata)
SE_testdata <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = mat_testdata))

test_that("incorrect input", {
  expect_error(
    runTBsigProfiler(list(mat_testdata),
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    paste0("Invalid input data type. Accepted input formats are",
           " matrix, data.frame, or SummarizedExperiment. ",
           "Your input: list")
  )
  expect_error(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "nothing", parallel.sz = 1),
    paste0("Invalid algorithm. Supported algorithms are GSVA, ",
           "ssGSEA, and ASSIGN")
  )
  expect_error(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("fakegene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    paste0("ERROR: all valid outputs are empty.")
  )
})

#Test matrix input
test_that("matrix input", {
  expect_error(
    runTBsigProfiler(mat_testdata, useAssay = "test",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "useAssay only supported for SummarizedExperiment objects"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ssGSEA", parallel.sz = 1),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata, outputFormat = "data.frame",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "data.frame"
  )
  expect_s4_class(
    runTBsigProfiler(mat_testdata, outputFormat = "SummarizedExperiment",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                     ASSIGNiter = 100, ASSIGNburnin = 50, parallel.sz = 1),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = FALSE,
                     algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                     ASSIGNiter = 100, ASSIGNburnin = 50, parallel.sz = 1),
    "matrix"
  )
})

#Test data.frame input
test_that("data.frame input", {
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ssGSEA", parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata, outputFormat = "matrix",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "matrix"
  )
  expect_s4_class(
    runTBsigProfiler(df_testdata, outputFormat = "SummarizedExperiment",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_is(
    runTBsigProfiler(df_testdata, outputFormat = "matrix",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = FALSE,
                     algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1),
    "data.frame"
  )
})

#Test SummarizedExperiment input
test_that("SummarizedExperiment input", {
  expect_error(
    runTBsigProfiler(SE_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "useAssay required for SummarizedExperiment Input"
  )
  expect_s4_class(
    runTBsigProfiler(SE_testdata, useAssay = "data",
                     algorithm = "GSVA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SE_testdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ssGSEA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SE_testdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SE_testdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = TRUE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment"
  )
  expect_error(
    runTBsigProfiler(SE_testdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = FALSE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment not supported with combineSigAndAlgorithm FALSE."
  )
})
