## code to prepare `DATASET` dataset goes here

usethis::use_data("DATASET")

library(purrr)
library(devtools)
############ Function to add a list of summarized experiment data ############
add_new_data <- function(edit_files){
  library(purrr)
  library(devtools)
  sobject_total <- lapply(edit_files,
                          function(x) readRDS(paste0("data-raw/",x)))
  sobject_total_name_ori <- lapply(sobject_total, function(x) names(assays(x)))
  sobject_total_name <- lapply(sobject_total_name_ori, function(x) gsub("(.+?)(\\_.*)", "\\1", x)) %>% unlist()
  sobject_total_name
  names(sobject_total) <- sobject_total_name

  purrr::walk2(sobject_total, paste0(names(sobject_total), "_sobject"), function(obj, name) {
    assign(name, obj)
    do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
  })

}

#----------------------------------------------------
# Add 3 data on Jan. 24th
edit_files <- paste0(c("GSE62147_sobject","GSE25534_sobject","GSE41055_sobject"),".RDS")
add_new_data(edit_files)

#----------------------------------------------------
# Modify RNA-seq data on Feb. 22th
add_new_data_mobject <- function(edit_files){
  library(purrr)
  library(devtools)
  sobject_total <- lapply(edit_files,
                          function(x) readRDS(paste0("~/Desktop/curatedTBData/data-raw/",x)))
  # sobject_total_name_ori <- lapply(sobject_total, function(x) names(assays(x)))
  # sobject_total_name <- lapply(sobject_total_name_ori, function(x) gsub("(.+?)(\\_.*)", "\\1", x)) %>% unlist()
  # sobject_total_name
  names(sobject_total) <- gsub("\\..*","",edit_files)

  purrr::walk2(sobject_total, paste0(names(sobject_total)), function(obj, name) {
    assign(name, obj)
    do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
  })

}

edit_files <- paste0(c("GSE101705_mobject", "GSE107104_mobject", "GSE89403_mobject",
                      "GSE79362_mobject","GSE94438_mobject", "GSE107991_mobject",
                      "GSE107992_mobject","GSE107993_mobject","GSE107994_mobject"),".RDS")
add_new_data_mobject(edit_files)

#----------------------------------------------------.
# Add 1 data on Mar. 1st
GSEXXXXX_sobject <- readRDS("data-raw/GSEXXXXX_sobject.RDS")
usethis::use_data(GSEXXXXX_sobject)

#----------------------------------------------------.
# Edit GSE112104 on Apr. 1st
add_new_data_mobject("GSE112104_mobject.RDS")

#----------------------------------------------------.
# Edit GSE107994,GSE79362 on Apr. 5th
add_new_data_mobject("GSE107994_mobject.RDS")
add_new_data_mobject("GSE79362_mobject.RDS")

#----------------------------------------------------.
# Edit GSE94438 on Apr. 5th
add_new_data_mobject("GSE94438_mobject.RDS")

#----------------------------------------------------.
# store Objects into indivual matrix

f1 <-  list.files("data-raw")
f2 <-  f1[grep("GSE",f1)]
total <- lapply(f2, function(x) readRDS(paste0("data-raw/",x)))

names(total) <- gsub("\\..*","",f2)

library(devtools)
purrr::walk2(total, names(total), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

# Save Summary table
DataSummary <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                  sheet = "DataSummary")
use_data(DataSummary,compress = "xz", overwrite = TRUE)

SignatureInfo <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                  sheet = "SignatureInfo")
use_data(SignatureInfo,compress = "xz", overwrite = TRUE)

SignatureInfoTraining <- readxl::read_excel("data-raw/Data_summaryforpackage.xlsx",
                                    sheet = "SignatureInfoTraining")
use_data(SignatureInfoTraining,compress = "xz", overwrite = TRUE)

library(devtools)
GSE79362_assay_reprocess <- GSE79362_raw_counts1
use_data(GSE79362_assay_reprocess,compress = "xz", overwrite = TRUE)
saveRDS(GSE79362_assay_reprocess,"/data-raw/GSE79362_assay_reprocess.RDS")


# Example to save a list of files from the same GEO
GSE22098_files <- list.files("data-raw",pattern = "GSE22098")
GSE22098_list <- lapply(GSE22098_files, function(x)
  readRDS(paste0("data-raw/",x)))

names(GSE22098_list) <- gsub("\\..*","",GSE22098_files)
library(devtools)
purrr::walk2(GSE22098_list, names(GSE22098_list), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

save_files_list <- function(geo_access){
  GSE_files <- list.files("data-raw",pattern = geo_access)
  GSE_list <- lapply(GSE_files, function(x)
    readRDS(paste0("data-raw/",x)))

  names(GSE_list) <- gsub("\\..*","",GSE_files)
  library(devtools)
  purrr::walk2(GSE_list, names(GSE_list), function(obj, name) {
    assign(name, obj)
    do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
  })
}

save_files_list("GSEIndia")
save_files_list("GSE83456")
save_files_list("GSE81746")
save_files_list("GSE107731")
save_files_list("GSE29536")

save_files_list("GSE83892")
get_curatedTBData("GSE84076")

# file_name <- c("GSE107991", "GSE107992", "GSE107993", "GSE107994")
file_name <- c("GSE94438")
GSE_list <- lapply(file_name, function(x)
  readRDS(paste0("data-raw/",x,"_assay_reprocess.RDS")))

names(GSE_list) <- paste0(file_name,"_assay_reprocess")
library(devtools)
purrr::walk2(GSE_list, names(GSE_list), function(obj, name) {
  assign(name, obj)
  do.call("use_data", list(as.name(name), compress = "xz", overwrite = TRUE))
})

# object_sub <- get_curatedTBData(geo_access = file_name)
# object_norm_test <- bplapply(object_sub, function(x)
#   Normalization(x, microarray_method = "quantile", RNAseq_method = "TMM",
#                 experiment_name = "assay_raw"), BPPARAM = param)
# object_match_test <- bplapply(object_norm_test, function(x)
#   MatchProbe(x, useAssay = c("TMM","quantile","RMA"),
#              createExperimentName = "assay_MatchProbe"), BPPARAM = param)

