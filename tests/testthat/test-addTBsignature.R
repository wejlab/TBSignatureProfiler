context("test-addTBsignature")


test_that("Run addTBsignature", {
  expect_message(
    addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                   authname = "Katarvi", signame_common = NULL,
                   sigtype = "Disease", tissuetype = "PBMC",
                   saveobjs = FALSE, views = FALSE)
  )
  expect_message(
    addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                   authname = "Katarvi", signame_common = "KAT3",
                   sigtype = "Disease/HIV", tissuetype = "whole blood",
                   saveobjs = FALSE, views = FALSE)
  )
})

test_that("Warnings", {
  expect_error(addTBsignature(sigsymbols = c(1, 2, 3),
                              authname = "Katarvi", signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("'sigsymbols' must be a character vector.")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                              authname = 1, signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("'authname' must be a character string.")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1"),
                              authname = "Katarvi", signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("Currently the TBSP only accepts gene signatures with 2",
                     "or more gene symbols.")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                              authname = "Katarvi", signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "PBMC",
                              saveobjs = NULL, views = FALSE),
               paste("'saveobjs' must be logical.")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                              authname = "Katarvi", signame_common = "KAT3",
                              sigtype = "Response", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("Input string for 'sigtype' not recognized.",
                     "You must denote 'sigtype' as one of the following:",
                     "'Disease', 'Disease/HIV', 'Disease/Other Diseases',",
                     "'Disease/Pneumonia', 'failure', 'response', 'risk'")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                              authname = "Katarvi", signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "Not",
                              saveobjs = FALSE, views = FALSE),
               paste("'tissuetype' must be one of the following:",
                     "'mixed', 'PBMC', 'whole blood'")
  )
  expect_error(addTBsignature(sigsymbols = c("APOL1", "BATF2", "ETV7"),
                              authname = "Thisisareallylongnamesok",
                              signame_common = "KAT3",
                              sigtype = "Disease/HIV", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("Resulting signature name is too long.",
                     "Please shorten or abbreviate authname to", 13,
                     "characters.")
  )
  expect_error(addTBsignature(sigsymbols = rep(c("APOL1", "BATF2"), 8),
                              authname = "Zak",
                              signame_common = "KAT3",
                              sigtype = "risk", tissuetype = "PBMC",
                              saveobjs = FALSE, views = FALSE),
               paste("There is already a signature with that name.",
                     "Please alter the authname to distinguish it from the previous",
                     "signature.")
  )
})
