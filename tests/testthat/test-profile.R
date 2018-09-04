context("test-profile")

#Test matrix input
test_that("matrix input", {
  expect_equal(2 * 2, 4)
})

#Test data.frame input
test_that("data.frame input", {
  expect_equal(2 * 2, 4)
})

#Test SummarizedExperiment input
test_that("SummarizedExperiment input", {
  expect_equal(2 * 2, 4)
})

#Test SingleCellExperiment input
test_that("SingleCellExperiment input", {
  expect_equal(2 * 2, 4)
})

#Test SCtkExperiment input
test_that("SCtkExperiment input", {
  expect_equal(2 * 2, 4)
})
