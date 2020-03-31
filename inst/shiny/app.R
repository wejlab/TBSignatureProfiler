if("shiny" %in% rownames(installed.packages()) == FALSE) {install.packages("shiny")}
if("shinyWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("shinyWidgets")}
if("shinythemes" %in% rownames(installed.packages()) == FALSE) {install.packages("shinythemes")}
if("markdown" %in% rownames(installed.packages()) == FALSE) {install.packages("markdown")}
if("DT" %in% rownames(installed.packages()) == FALSE) {install.packages("DT")}
if("BiocManager" %in% rownames(installed.packages()) == FALSE) {install.packages("BiocManager")}
if("SummarizedExperiment" %in% rownames(installed.packages()) == FALSE) {install.packages("SummarizedExperiment")}
if("TBSignatureProfiler" %in% rownames(installed.packages()) == FALSE) {devtools::install_github("compbiomed/TBSignatureProfiler")}
if("shinycssloaders" %in% rownames(installed.packages()) == FALSE) {install.packages("shinycssloaders")}

library(shiny)
library(shinyWidgets)
library(shinythemes)
library(markdown)
library(DT)
library(BiocManager)
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(shinycssloaders)


siglist=names(TBsignatures)

ui <- fluidPage(theme=shinytheme('cosmo'),
                navbarPage("TB Signature Profiler",
                           source(file.path('ui', 'ui_upload.R'), local=T)$value,
                           source(file.path('ui', 'ui_profiler.R'), local=T)$value,
                           source(file.path('ui', 'ui_visualization.R'), local=T)$value,
                           source(file.path('ui', 'ui_auc.R'), local=T)$value
                )


)

server <- function(input, output, session) {

  vals=reactiveValues(
    tbdat=NULL,
    coldat=NULL,
    covars=NULL,
    datassays=NULL,
    profilerdat=NULL
  )

  source(file.path('server', 'server_upload.R'), local=T)$value
  source(file.path('server', 'server_profiler.R'), local=T)$value
  source(file.path('server', 'server_visualization.R'), local=T)$value
  source(file.path('server', 'server_auc.R'), local=T)$value

}

shinyApp(ui, server)
