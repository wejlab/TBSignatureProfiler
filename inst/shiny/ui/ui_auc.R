shiny::navbarMenu(
  "AUCs Via Bootstrapping",
  shiny::tabPanel(
    "T-tests & AUC",
    DT::dataTableOutput("ttest") %>% shinycssloaders::withSpinner()
  ),
  shiny::tabPanel(
    "Boxplots for the bootstrapped AUCs",
    shiny::sidebarPanel(
      shinyWidgets::pickerInput(
        "bootsigs", "Select Signatures", choices = siglist,
        options = list("actions-box" = TRUE), multiple = TRUE, selected = NULL),
      shiny::actionButton("bootplot", "Plot Bootstrapped AUCs")
    ),
    shiny::mainPanel(plotOutput("bootbox", height = 600))
  ),
  shiny::tabPanel(
    "Separate ROC plots, 95% CI Bands",
    shiny::sidebarPanel(
      shinyWidgets::pickerInput(
        "singroc", "Select Signature(s)", choices = siglist,
        options = list("actions-box" = TRUE), multiple = TRUE, selected = NULL),
      shiny::actionButton("rocplot", "Plot ROC Curve")
    ),
    shiny::mainPanel(
      shiny::plotOutput("rocsep", height = 450)
    )
  )
)
