shiny::tabPanel(
  "Upload Data",
  tags$div(
    class = "jumbotron", style = "background-color:#ededed",
    tags$div(
      class = "container",
      shiny::h1("TB Signature Profiler"),
      shiny::p("Score, compare and visualize gene signatures"),
    )
  ),
  shiny::sidebarPanel(
    shiny::radioButtons(
      "dat", "Select Data", choices = c("Sample Data" = "samp",
                                        "Upload Data" = "updat")),
    shiny::conditionalPanel(
      condition = "input.dat == 'samp'",
      shiny::radioButtons(
        "sampdat", "Sample Datasets",
        choices = c("TB HIV" = "tbhiv", "TB Indian" = "tbind"),
        selected = NULL)
    ),
    shiny::conditionalPanel(
      condition = "input.dat == 'updat'",
      shiny::radioButtons(
        "updattype", "Choose File Type",
        choices = c(
          "Tab Separated Text Files" = "rawdat",
          "Summarized Experiment" = "sumexp"),
        selected = NULL)
    ),
    shiny::hr(),
    shiny::actionButton("upload", "Upload")),
  shiny::mainPanel(
    shiny::conditionalPanel(
      condition = "input.dat == 'updat'",
      shiny::conditionalPanel(
        condition = "input.updattype == 'rawdat'",
        shiny::column(
          width = 4,
          shiny::wellPanel(
            shiny::h4("Upload Count File"),
            shiny::fileInput(
              "countsfile",
              shiny::HTML(
                base::paste(
                  "Input assay (eg. counts, required):",
                  tags$span(style = "color:red", "*", sep = ""))
              ),
              accept = c(
                "text/csv",
                "text/comma-separated-values",
                "text/tab-separated-values",
                "text/plain", ".csv", ".tsv"
              )
            )
          )
        ),
        shiny::column(
          width = 4,
          shiny::wellPanel(
            shiny::h4("Upload Meta Data File"),
            shiny::fileInput(
              "metdatfile",
              shiny::HTML(
                base::paste(
                  "Meta Data (eg. disease, required):",
                  tags$span(style = "color:red", "*", sep = ""))
              ),
              accept = c(
                "text/csv",
                "text/comma-separated-values",
                "text/tab-separated-values",
                "text/plain", ".csv", ".tsv"
              )
            )
          )
        ),
        shiny::selectInput(
          "datassay",
          label = "Input Assay Type:",
          choices = c("Counts" = "counts", "Log Counts" = "log",
                      "CPM" = "cpm", "Log CPM" = "logcpm")
        )
      ),
      shiny::conditionalPanel(
        condition = "input.updattype == 'sumexp'",
        shiny::h3(paste("Choose an RDS Summarized Experiment file",
                        "that contains a TBSignature Object:")),
        shiny::fileInput("rdsFile", "TBSignature RDS file:",
                         accept = c(".rds", ".RDS"))
      )
    ),
    DT::dataTableOutput("counttable", width = 700)
  )
)
