#informal testing
library(TBSignatureProfiler)

india_relapse <- readRDS("~/Dropbox/grad_school/johnson_lab/20180816_failure_relapse_india_2018/failure_relapse_india_2018.rds")
india_relapse_GSVA <- runTBsigProfiler(input=india_relapse, useAssay = "logtpm", algorithm = c("GSVA", "ssGSEA", "ASSIGN"), assignDir = "~/Dropbox/grad_school/johnson_lab/20180816_failure_relapse_india_2018/assignRes", combineSigAndAlgorithm = TRUE)
india_relapse_GSVA_df <- runTBsigProfiler(input=india_relapse, useAssay = "logtpm", algorithm = c("GSVA", "ssGSEA", "ASSIGN"), assignDir = "~/Dropbox/grad_school/johnson_lab/20180816_failure_relapse_india_2018/assignRes", outputFormat = "data.frame", combineSigAndAlgorithm = TRUE)

#SHOULD NOT ERROR
# india_relapse_GSVA_df <- runTBsigProfiler(input=india_relapse, useAssay = "logtpm", algorithm = c("GSVA", "ssGSEA", "ASSIGN"), outputFormat = "data.frame", combineSigAndAlgorithm = F)
# Error in runTBsigProfiler(input = india_relapse, useAssay = "logtpm",  :
#                             SummarizedExperiment not supported with combineSigAndAlgorithm FALSE.

colData(india_relapse_GSVA)$visit_subject <- paste(colData(india_relapse_GSVA)$visit, colData(india_relapse_GSVA)$subjtype, sep="\n")

a <- signatureBoxplot(
  india_relapse_GSVA,
  signatureColNames = c("GSVA_ACS_COR", "ssGSEA_ACS_COR", "ASSIGN_ACS_COR"),
  annotationColName = c("visit_subject"),
  name = "ACS_COR Signatures")
a

b <- signatureHeatmap(
  india_relapse_GSVA,
  signatureColNames = colnames(colData(india_relapse_GSVA))[71:115],
  annotationColName = c("visit", "subjtype"), scale=TRUE,
  showColumnNames = FALSE)

b <- signatureGeneHeatmap(
  india_relapse_GSVA,
  useAssay = "logtpm",
  sigGenes = TBsignatures[["ACS_COR"]],
  signatureColNames = c("GSVA_ACS_COR", "ssGSEA_ACS_COR", "ASSIGN_ACS_COR"),
  annotationColNames = c("visit", "subjtype"),
  showColumnNames = FALSE,
  name = "ACS_COR")

#' #generate some artificial data that shows a difference in Predict29
#' mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
#'                              dimnames = list(TBsignatures$Predict29,
#'                                              paste0("sample", 1:10))),
#'                       matrix(rnorm(1000), 100, 10,
#'                              dimnames = list(paste0("gene", 1:100),
#'                                              paste0("sample", 1:10))))
#' #Create a SummarizedExperiment object that contains the data
#' testdataSE <- SummarizedExperiment(assays=SimpleList(data = mat_testdata),
#'                                    colData=DataFrame(sample = c(rep("down", 5),
#'                                                                 rep("up", 5))))
#' #Run profiler using GSVA and ssGSEA on Predict29
#' res <- runTBsigProfiler(testdataSE, useAssay = "data",
#'                         signatures = TBsignatures["Predict29"],
#'                         algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1,
#'                         combineSigAndAlgorithm = TRUE)
#' #Plot heatmap of the results
#' signatureHeatmap(res, signatureColNames = c("GSVA_Predict29",
#'                                             "ssGSEA_Predict29"),
#'                  annotationColName = "sample", scale=TRUE,
#'                  showColumnNames = FALSE)

#generate some artificial data that shows a difference in Predict29
mat_testdata <- rbind(matrix(c(rnorm(145), rnorm(145) + 5), 29, 10,
                             dimnames = list(TBsignatures$Predict29,
                                             paste0("sample", 1:10))),
                      matrix(rnorm(1000), 100, 10,
                             dimnames = list(paste0("gene", 1:100),
                                             paste0("sample", 1:10))))
#Create a SummarizedExperiment object that contains the data
testdataSE <- SummarizedExperiment(assays=SimpleList(data = mat_testdata),
                                   colData=DataFrame(sample = c(rep("down", 5),
                                                                rep("up", 5))))
#Run profiler using GSVA and ssGSEA on Predict29
res <- runTBsigProfiler(testdataSE, useAssay = "data",
                        signatures = TBsignatures["Predict29"],
                        algorithm = c("GSVA", "ssGSEA"), parallel.sz = 1,
                        combineSigAndAlgorithm = TRUE)
#Plot heatmap of the results
signatureHeatmap(res, signatureColNames = c("GSVA_Predict29",
                                            "ssGSEA_Predict29"),
                 annotationColNames = "sample", scale=TRUE,
                 showColumnNames = FALSE)

#' signatureBoxplot(res, signatureColNames = c("GSVA_Predict29",
#'                                             "ssGSEA_Predict29"),
#'                  annotationColName = "sample", name = "Predict29 Signatures")

#' signatureGeneHeatmap(res, useAssay = "data",
#'                      sigGenes = TBsignatures[["Predict29"]],
#'                      signatureColNames = c("GSVA_Predict29",
#'                                            "ssGSEA_Predict29"),
#'                      annotationColNames = c("sample"), showColumnNames = FALSE,
#'                      name = "Predict29")

#TODO: input1: SummarizedExperiment object, list of columns of signature data, list of columns of annotation data
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 signatureColNames = c("GSVA_ACS_COR", "GSVA_Predict29", "ssGSEA_ACS_COR", "ssGSEA_Predict29", "ASSIGN_ACS_COR", "ASSIGN_Predict29"),
                 annotationColNames = c("visit", "subjtype"),
                 scale=FALSE)
#TODO: input2: signature data, annotation data frame (matrix of signature data and matrix of annotations)
annottest <- data.frame(colData(india_relapse_GSVA)[,c("visit", "subjtype")])
signatureHeatmap(inputData = india_relapse_GSVA_df,
                 annotationData = annottest,
                 name="TEST SIGS",
                 scale=FALSE)
#TODO: try with one signature
#TODO: try with one annotation
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 signatureColNames = c("GSVA_ACS_COR"),
                 annotationColNames = c("subjtype"),
                 scale=FALSE)
signatureHeatmap(inputData = india_relapse_GSVA_df[12,,drop=F],
                 annotationData = annottest[,2,drop=F],
                 name="TEST SIGS",
                 scale=FALSE)
#TODO: check if you have a duplicate name in colnames values
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 signatureColNames = c("GSVA_ACS_COR","GSVA_ACS_COR","GSVA_ACS_COR"),
                 annotationColNames = c("subjtype"),
                 scale=FALSE)
# signatureHeatmap(india_relapse_GSVA, signatureColNames = colnames(colData(india_relapse_GSVA))[71:85],
# annotationColNames=c("visit", "subjtype"))

#check if you don't provide annotationColNames
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 signatureColNames = c("GSVA_ACS_COR"),
                 scale=FALSE)
#check if you don't provide sigColNames
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 scale=FALSE)
#duplicate signature column names
signatureHeatmap(inputData = india_relapse_GSVA,
                 annotationData = FALSE,
                 name="TEST SIGS2",
                 signatureColNames = c("GSVA_ACS_COR","GSVA_ACS_COR","GSVA_ACS_COR"),
                 scale=FALSE)


