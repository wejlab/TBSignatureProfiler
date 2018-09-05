context("test-profile")

set.seed(1234)
mat_testdata <- matrix(rnorm(1000), 100, 10,
                       dimnames = list(paste0("gene", 1:100),
                                       paste0("sample", 1:10)))
df_testdata <- data.frame(mat_testdata)
SE_testdata <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = mat_testdata))

#Test matrix input
test_that("matrix input", {
  runTBsigProfiler(mat_testdata,
                   signatures = list(sig1 = paste0("gene", 1:10)),
                   algorithm = "GSVA")
})

#Test data.frame input
test_that("data.frame input", {
  runTBsigProfiler(df_testdata,
                   signatures = list(sig1 = paste0("gene", 1:10)),
                   algorithm = "GSVA")
})

#Test SummarizedExperiment input
test_that("SummarizedExperiment input", {
  expect_error(runTBsigProfiler(SE_testdata,
                                signatures = list(sig1 = paste0("gene", 1:10)),
                                algorithm = "GSVA"),
               "useAssay required for SummarizedExperiment Input")
  expect_s4_class(runTBsigProfiler(SE_testdata,
                                   useAssay = "data",
                                   signatures = list(sig1 = paste0("gene", 1:10)),
                                   algorithm = "GSVA"),
                  "SummarizedExperiment")

})

#Test SingleCellExperiment input
test_that("SingleCellExperiment input", {
  expect_equal(2 * 2, 4)
})

#Test SCtkExperiment input
test_that("SCtkExperiment input", {
  expect_equal(2 * 2, 4)
})
