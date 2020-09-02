for (package in c("shiny", "shinyWidgets", "shinythemes", "markdown", "DT",
                  "BiocManager", "SummarizedExperiment", "shinycssloaders")) {
  if (identical(find.package(package, quiet = TRUE), character(0))) {
    install.packages(package)
  }
}

library(shiny)
library(shinyWidgets)
library(shinythemes)
library(markdown)
library(DT)
library(BiocManager)
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(shinycssloaders)

siglist <- names(TBsignatures)

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  navbarPage(
    "TB Signature Profiler",
    source(file.path("ui", "ui_upload.R"), local = TRUE)$value,
    source(file.path("ui", "ui_profiler.R"), local = TRUE)$value,
    source(file.path("ui", "ui_visualization.R"), local = TRUE)$value,
    source(file.path("ui", "ui_auc.R"), local = TRUE)$value
  )
)

server <- function(input, output, session) {
  vals <- reactiveValues(
    tbdat = NULL,
    coldat = NULL,
    covars = NULL,
    datassays = NULL,
    profilerdat = NULL
  )
  source(file.path("server", "server_upload.R"), local = TRUE)$value
  source(file.path("server", "server_profiler.R"), local = TRUE)$value
  source(file.path("server", "server_visualization.R"), local = TRUE)$value
  source(file.path("server", "server_auc.R"), local = TRUE)$value
}

shinyApp(ui, server)
