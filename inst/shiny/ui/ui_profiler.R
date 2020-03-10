tabPanel("Run TB Signature Profiler",
         sidebarPanel(
           pickerInput("newassay", label = "Make Assay",
                       choices = c("Log Counts"='log', "CPM"='cpm', "Log CPM"='logcpm')),
           actionButton("mkassay", "Make Assay"),
           hr(),
           selectInput("profassay", "Select Assay for Profiler", choices = NULL),
           selectInput('profalg', "Select Algorithm for Profiler", choices = c("GSVA", "ssGSEA", "ASSIGN", "PLAGE", "Zscore", "singscore")),
           actionButton('runprofiler', "Run Profiler")
         ),
         mainPanel(
           h4("Covariates"),
           tableOutput("covars"),
           h4("Assays"),
           tableOutput("assays"),
           h4("TBSignatureProfiler Data"),
           hr(),
           dataTableOutput("visdat", width = 700)
         )
)
