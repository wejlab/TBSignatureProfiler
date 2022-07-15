context("test-OriginalModel")
## Use "TB_hiv" as example for the input
mat_testdata <- rbind(matrix(c(rnorm(80), rnorm(80) + 5), 16, 10,
                             dimnames = list(TBsignatures$Zak_RISK_16,
                                             paste0("sample", 1:10))),
                      matrix(rnorm(1000), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
col_info <- S4Vectors::DataFrame(condition = rep(c("Disease", "Control"), each = 5),
                                 row.names = colnames(mat_testdata))
sobject <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat_testdata),
                                                      colData = col_info)

test_that("\"input\" should be a SummarizedExperiment object.", {
    expect_error(evaluateOriginalModel(input = mat_testdata,
                                       geneSignaturesName = "Zak_RISK_16"))
})

test_that("Argument \"geneSignaturesName\" cannot be missing.", {
    expect_error(evaluateOriginalModel(input = TB_hiv, useAssay = "counts"))
})

test_that("Gene sets not found from the \"input\"", {
    re <- suppressWarnings(evaluateOriginalModel(input = sobject,
                                geneSignaturesName = c("Anderson_42", "Leong_RISK_29"),
                                useAssay = "counts"))
    expect_identical(re$Leong_RISK_29_OriginalModel, rep(NA, dim(re)[2]))
    expect_identical(re$Anderson_42_OriginalModel, rep(NA, dim(re)[2]))
    expect_warning(evaluateOriginalModel(input = sobject,
                                       geneSignaturesName = "Leong_RISK_29",
                                       useAssay = "counts"))
    expect_warning(evaluateOriginalModel(input = sobject,
                                         geneSignaturesName = "Anderson_42",
                                         useAssay = "counts"))
})

test_that("Input signatures' original model not available.", {
    expect_error(evaluateOriginalModel(input = TB_hiv,
                                       geneSignaturesName = "EXAMPLE_NAME",
                                       useAssay = "counts"))
    expect_message(evaluateOriginalModel(input = TB_hiv,
                                         geneSignaturesName = c("EXAMPLE_NAME",
                                                                "Zak_RISK_16"),
                                         useAssay = "counts"))
})

test_that("Assay not found from the input.", {
    expect_error(evaluateOriginalModel(input = TB_hiv,
                                       geneSignaturesName = "Zak_RISK_16",
                                       useAssay = "EXAMPLE_ASSAY"))
})

test_that("Return type is a SummarizedExperiment object.", {
    re <- evaluateOriginalModel(input = TB_hiv,
                                geneSignaturesName = "Zak_RISK_16",
                                useAssay = "counts")
    expect_s4_class(re, "SummarizedExperiment")
})
