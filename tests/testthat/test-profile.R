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
    paste0("Invalid algorithm. Supported algorithms are: ",
           "GSVA, ssGSEA, PLAGE, Zscore, singscore, and ASSIGN")
  )
  expect_error(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("fakegene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    paste0("ERROR: all valid outputs are empty.")
  )
})

# Test matrix input
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
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "PLAGE"),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "Zscore"),
    "matrix"
  )
  expect_is(
    runTBsigProfiler(mat_testdata,
                     signatures = TBsignatures["Zak_RISK_16"],
                     algorithm = "singscore"),
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
  expect_error(
    runTBsigProfiler(mat_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = NULL,
                     algorithm = c("GSVA", "ssGSEA", "ASSIGN"),
                     ASSIGNiter = 100, ASSIGNburnin = 50, parallel.sz = 1),
    paste0("You must choose whether or not to combine the ",
           "signature and algorithm name using combineSigAndAlgorithm.")
  )
})

# Test data.frame input
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
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "PLAGE", ASSIGNiter = 100, ASSIGNburnin = 50),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "Zscore", ASSIGNiter = 100, ASSIGNburnin = 50),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = TBsignatures["Zak_RISK_16"],
                     algorithm = "singscore", ASSIGNiter = 100, ASSIGNburnin = 50),
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
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = FALSE,
                     algorithm = c("ssGSEA", "ASSIGN"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("ssGSEA", "ASSIGN"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("GSVA", "Zscore"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("ssGSEA", "Zscore"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("ASSIGN", "Zscore"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("GSVA", "PLAGE"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("ssGSEA", "PLAGE"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("ASSIGN", "PLAGE"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  expect_is(
    runTBsigProfiler(df_testdata,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     combineSigAndAlgorithm = TRUE,
                     algorithm = c("Zscore", "PLAGE"), ASSIGNiter = 100,
                     ASSIGNburnin = 50, parallel.sz = 1),
    "data.frame"
  )
  assignDir <- tempfile("assign")
  dir.create(assignDir)
  expect_is(
    runTBsigProfiler(df_testdata, assignDir = assignDir,
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "ASSIGN", ASSIGNiter = 100, ASSIGNburnin = 50),
    "data.frame"
  )
  unlink(assignDir)
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
                     signatures = TBsignatures["Zak_RISK_16"],
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
                     algorithm = "PLAGE", ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = "Zscore", ASSIGNiter = 100, ASSIGNburnin = 50),
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
    "SummarizedExperiment not supported when combineSigAndAlgorithm FALSE."
  )
  expect_error(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     outputFormat = "SummarizedExperiment",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA", "ASSIGN"), parallel.sz = 1,
                     combineSigAndAlgorithm = FALSE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "SummarizedExperiment not supported when combineSigAndAlgorithm FALSE."
  )
})

# Test other things
SEtestdata <- SummarizedExperiment::SummarizedExperiment(
  assays = S4Vectors::SimpleList(data = mat_testdata,
                                 counts = mat_testdata))
test_that("Output error", {
  expect_error(
    runTBsigProfiler(SEtestdata, useAssay = "data",
                     outputFormat = "other",
                     signatures = list(sig1 = paste0("gene", 1:10)),
                     algorithm = c("GSVA"), parallel.sz = 1,
                     combineSigAndAlgorithm = TRUE,
                     ASSIGNiter = 100, ASSIGNburnin = 50),
    "Output format error."
  )
  expect_is(
    update_genenames(c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4", "4-Oct",
                       "OCT4-PG4", "C19ORF71", "C19orf71")),
    "character"
  )
  expect_s4_class(
    runTBsigProfiler(SEtestdata, useAssay = NULL,
                     signatures = TBsignatures["Zak_RISK_16"],
                     algorithm = "GSVA", parallel.sz = 1),
    "SummarizedExperiment"
  )
})

# Test compareAlgs function
SEtestdata$Disease <- factor(c(rep("yes", 5), rep("no", 5)))
test_that("CompareAlgs function", {
  expect_error(
    compareAlgs(SEtestdata,
                useAssay = "data",
                annotationColName = "Disease",
                algorithm = "GSVA",
                output = "None",
                signatures = list(sig1 = paste0("gene", 1:10))),
    "Output parameter must specify either 'heatmap' or 'boxplot'"
  )
  expect_is(
    compareAlgs(SEtestdata,
                useAssay = "data",
                annotationColName = "Disease",
                algorithm = "GSVA",
                output = "heatmap",
                signatures = list(sig1 = paste0("gene", 1:10))),
    "HeatmapList"
  )
  expect_output(
    compareAlgs(SEtestdata,
                useAssay = "data",
                annotationColName = "Disease",
                algorithm = "GSVA",
                output = "boxplot",
                show.pb = TRUE,
                signatures = list(sig1 = paste0("gene", 1:5),
                                  sig2 = paste0("gene", 5:10)))
  )
  expect_error(
    compareAlgs(df_testdata,
                useAssay = NULL,
                annotationColName = "Disease",
                algorithm = "GSVA",
                output = "boxplot",
                show.pb = TRUE,
                signatures = list(sig1 = paste0("gene", 1:10))),
    "Input must be a SummarizedExperiment object."
  )
})
