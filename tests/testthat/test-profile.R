context("test-profile")

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
    paste0("Invalid algorithm. Supported algorithms are
    GSVA, ssGSEA, PLAGE, Zscore, singscore, and ASSIGN")
  )
  expect_error(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("fakegene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    paste0("ERROR: all valid outputs are empty.")
  )
})

#Test SummarizedExperiment input
test_that("SummarizedExperiment input", {
  expect_error(
    runTBsigProfiler(SEtestdata, signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "GSVA", parallel.sz = 1),
    "useAssay required for SummarizedExperiment Input"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = TBsignatures['Zak_RISK_16'],
                     algorithm = "GSVA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ssGSEA", parallel.sz = 1),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = TRUE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment"
  )
  expect_error(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = FALSE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment not supported with combineSigAndAlgorithm FALSE."
  )
  expect_error(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     outputFormat = "SummarizedExperiment",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = FALSE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment not supported with combineSigAndAlgorithm FALSE."
  )
})
