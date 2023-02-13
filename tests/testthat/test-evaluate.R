context("test-evaluate")

num_sample <- 20
num_genes <- 100
inputTestmat <- matrix(rnorm(1000), num_genes, num_sample,
                    dimnames = list(paste0("gene", seq.int(1, num_genes)),
                                    paste0("sample", seq.int(1, num_sample))))
inputTest <- as.data.frame(inputTestmat)
target_Vec <- sample(c(0, 1), replace = TRUE, size = num_sample)

PQ_target_Vec <- c(rep(1, num_sample/2), rep(0, num_sample/2))
signaturelist <- list(sig1 = c("gene1", "gene2", "gene3"),
                       sig2 = c("gene4", "gene5", "gene6"))
signaturenamevec <- c("sig1", "sig2")
numboot <- 5

test_that("deseq2_norm_rle works", {
  expect_is(
    deseq2_norm_rle(abs(inputTest)),
    "data.frame"
  )
})

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
                        pb.show = FALSE),
  "list"
  )
  expect_error(
    SignatureQuantitative(inputTest,
                          targetVec.num = target_Vec,
                          signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
                                                sig2 = c("gene4", "gene5", "gene6")),
                          signature.name.vec = NULL,
                          num.boot = 3,
                          pb.show = FALSE),
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
                          pb.show = FALSE),
    "The inputs signature.list and signature.name.vec are not the same
         length."
  )
})

#test_that("plotQuantitative works", {
#  expect_is(
#    plotQuantitative(inputTest,
#                     targetVec.num = PQ_target_Vec,
#                     signature.list = signaturelist,
#                     signature.name.vec = signaturenamevec,
#                     num.boot = 5,
#                     pb.show = FALSE,
#                     rotateLabels = TRUE),
#    "gg"
#  )
#  expect_error(
#    plotQuantitative(inputTest,
#                     targetVec.num = PQ_target_Vec,
#                     signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
#                                           sig2 = c("gene4", "gene5", "gene6")),
#                     signature.name.vec = NULL,
#                     num.boot = 3,
#                     pb.show = FALSE,
#                     rotateLabels = TRUE),
#    "Please specify arguments for both signature.list and
#         signature.name.vec, or leave them both empty to use
#         TBsignatures as the list of signatures for profiling."
#  )
#  expect_error(
#    plotQuantitative(inputTest,
#                     targetVec.num = PQ_target_Vec,
#                     signature.list = list(sig1 = c("gene1", "gene2", "gene3"),
#                                           sig2 = c("gene4", "gene5", "gene6")),
#                     signature.name.vec = c("sig2"),
#                     num.boot = 3,
#                     pb.show = FALSE,
#                     rotateLabels = TRUE),
#    "The inputs signature.list and signature.name.vec are not the same
#         length."
#  )
#})
#