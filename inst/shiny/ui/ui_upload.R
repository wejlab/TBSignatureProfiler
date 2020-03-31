tabPanel("Upload Data",
         tags$div(
           class = "jumbotron", style = "background-color:#ededed",
           tags$div(
             class = "container",
             h1("TB Signature Profiler"),
             p("Score, compare and visualize gene signatures"),
           )
         ),
         sidebarPanel(
           radioButtons("dat", "Select Data", choices = c("Sample Data"="samp", "Upload Data"="updat")),
           conditionalPanel(condition = "input.dat=='samp'",
                            radioButtons("sampdat", "Sample Datasets", choices = c('TB HIV'='tbhiv', 'TB Indian'='tbind'), selected = NULL)
           ),
           conditionalPanel(condition = "input.dat=='updat'",
                            radioButtons("updattype","Choose File Type", choices=c('Tab Separated Text Files'='rawdat', 'Summarized Experiment'='sumexp'), selected = NULL)
           ),
           hr(),
           actionButton('upload', 'Upload')
         ),
         mainPanel(
           conditionalPanel(condition="input.dat=='updat'",
                            conditionalPanel(condition= "input.updattype=='rawdat'",
                                             column(width = 4,
                                                    wellPanel(
                                                      h4("Upload Count File"),
                                                      fileInput(
                                                        "countsfile",
                                                        HTML(
                                                          paste("Input assay (eg. counts, required):",
                                                                tags$span(style = "color:red", "*", sep = ""))
                                                        ),
                                                        accept = c(
                                                          "text/csv", "text/comma-separated-values",
                                                          "text/tab-separated-values", "text/plain", ".csv", ".tsv"
                                                        )
                                                      )
                                                    )
                                             ),
                                             column(width = 4,
                                                    wellPanel(
                                                      h4("Upload Meta Data File"),
                                                      fileInput(
                                                        "metdatfile",
                                                        HTML(
                                                          paste("Meta Data (eg. disease, required):",
                                                                tags$span(style = "color:red", "*", sep = ""))
                                                        ),
                                                        accept = c(
                                                          "text/csv", "text/comma-separated-values",
                                                          "text/tab-separated-values", "text/plain", ".csv", ".tsv"
                                                        )
                                                      )
                                                    )
                                             ),
                                             selectInput("datassay", label = "Input Assay Type:",
                                                         choices = c("Counts"="counts","Log Counts"='log', "CPM"='cpm', "Log CPM"='logcpm'))
                            ),
                            conditionalPanel(condition= "input.updattype=='sumexp'",
                                             h3("Choose an RDS Summarized Experiment file that contains a TBSignature Object:"),
                                             fileInput(
                                               "rdsFile", "TBSignature RDS file:", accept = c(".rds", ".RDS")
                                             )
                            )
           ),
           DT::dataTableOutput("counttable", width = 700)
         )
)
