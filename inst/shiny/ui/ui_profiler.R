shiny::tabPanel("Run TB Signature Profiler",
                shiny::sidebarPanel(
                  shinyWidgets::pickerInput("newassay", label = "Make Assay",
                                            choices = c(
                                              "Log Counts" = 'log',
                                              "CPM" = 'cpm',
                                              "Log CPM" = 'logcpm')),
                  shiny::actionButton("mkassay", "Make Assay"),
                  shiny::hr(),
                  shiny::selectInput("profassay", "Select Assay for Profiler",
                                     choices = NULL),
                  shiny::selectInput('profalg', "Select Algorithm for Profiler",
                                     choices = c(
                                       "GSVA", "ssGSEA", "ASSIGN", "PLAGE",
                                       "Zscore", "singscore")),
                  shiny::actionButton('runprofiler', "Run Profiler")
                  ),
                shiny::mainPanel(
                  shiny::h4("Covariates"),
                  shiny::tableOutput("covars"),
                  shiny::h4("Assays"),
                  shiny::tableOutput("assays"),
                  shiny::h4("TBSignatureProfiler Data"),
                  shiny::hr(),
                  DT::dataTableOutput("visdat", width = 700)
                  )
                )
