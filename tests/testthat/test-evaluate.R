context("test-evaluate")

set.seed(1234)
inputTestmat <- matrix(rnorm(1000), 100, 20,
                    dimnames = list(paste0("gene", seq.int(1, 100)),
                                    paste0("sample", seq.int(1, 20))))
inputTest <- as.data.frame(inputTestmat)
target_Vec <- sample(c(0, 1), replace = TRUE, size = 20)
signaturelist <- list(sig1 = c("gene1", "gene2", "gene3"),
                       sig2 = c("gene4", "gene5", "gene6"))
signaturenamevec <- c("sig1", "sig2")
numboot <- 5

test_that("LOOAUC works", {
  expect_is(
    LOOAUC_simple_multiple_noplot_one_df(inputTest,
                                         targetVec = target_Vec),
    "list"
  )
  expect_is(
    Bootstrap_LOOCV_LR_AUC(inputTest,
                           targetVec = target_Vec,
                           nboot = 10),
    "list"
  )
})

test_that("signatureQuantitative works", {
  expect_is(
  SignatureQuantitative(inputTest,
                        targetVec.num = target_Vec,
                        signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                               sig2 = c("gene4", "gene5", "gene6")),
                        signature.name.vec = c("sig1", "sig2"),
                        num.boot = 3,
                        pb.show = TRUE),
  "list"
  )
  expect_error(
    SignatureQuantitative(inputTest,
                          targetVec.num = target_Vec,
                          signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                                sig2 = c("gene4", "gene5", "gene6")),
                          signature.name.vec = NULL,
                          num.boot = 3,
                          pb.show = TRUE),
    "Please specify arguments for both signature.list and
         signature.name.vec, or leave them both empty to use
         TBsignatures as the list of signatures for profiling."
  )
  expect_error(
    SignatureQuantitative(inputTest,
                          targetVec.num = target_Vec,
                          signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                                sig2 = c("gene4", "gene5", "gene6")),
                          signature.name.vec = c("sig1"),
                          num.boot = 3,
                          pb.show = TRUE),
    "The inputs signature.list and signature.name.vec are not the same
         length."
  )
})

test_that("plotQuantitative works", {
  expect_is(
    plotQuantitative(inputTest,
                     targetVec.num = target_Vec,
                     signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                           sig2 = c("gene4", "gene5", "gene6")),
                     signature.name.vec = c("sig1", "sig2"),
                     num.boot = 20,
                     pb.show = TRUE,
                     rotateLabels = TRUE),
    "gg"
  )
  expect_error(
    plotQuantitative(inputTest,
                     targetVec.num = target_Vec,
                     signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                           sig2 = c("gene4", "gene5", "gene6")),
                     signature.name.vec = NULL,
                     num.boot = 3,
                     pb.show = TRUE,
                     rotateLabels = TRUE),
    "Please specify arguments for both signature.list and
         signature.name.vec, or leave them both empty to use
         TBsignatures as the list of signatures for profiling."
  )
  expect_error(
    plotQuantitative(inputTest,
                     targetVec.num = target_Vec,
                     signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                           sig2 = c("gene4", "gene5", "gene6")),
                     signature.name.vec = c("sig2"),
                     num.boot = 3,
                     pb.show = TRUE,
                     rotateLabels = TRUE),
    "The inputs signature.list and signature.name.vec are not the same
         length."
  )
})
