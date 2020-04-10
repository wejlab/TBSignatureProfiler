shiny::navbarMenu("Signature Visualization",
                  shiny::tabPanel("Heatmap with Selected Signatures",
                                  shiny::sidebarPanel(
                                    shinyWidgets::pickerInput(
                                      'selectsigs',
                                      label = 'Select Signatures for Profiler',
                                      choices = siglist,
                                      options = list('actions-box' = TRUE),
                                      multiple = TRUE, selected = siglist),
                                    shiny::selectInput(
                                      'allheatcovar',
                                      "Select Covariate",
                                      choices = NULL),
                                    shiny::actionButton('allheatplot',
                                                        "Plot Heatmap")
                                    ),
                                  shiny::mainPanel(
                                    plotOutput("allheat", height = 800)
                                    )
                                  ),
                  shiny::tabPanel("Heatmaps of Individual Signatures",
                                  shiny::sidebarPanel(
                                    shinyWidgets::pickerInput(
                                      'singheat',
                                      'Signature',
                                      choices = siglist,
                                      multiple = FALSE,
                                      selected = NULL),
                                    shinyWidgets::pickerInput(
                                      'genes', label = 'Genes',
                                      choices = NULL,
                                      options = list(
                                        'actions-box' = TRUE),
                                      multiple = TRUE),
                                    shiny::hr(),
                                    shiny::selectInput('singheatcovar',
                                                       'Covariate',
                                                       choices = NULL),
                                    shiny::actionButton(
                                      'singheatplot', "Plot Heatmap(s)")
                                    ),
                                  shiny::mainPanel(
                                    shiny::plotOutput("indheat", height = 800)
                                    )
                                  ),
                  shiny::tabPanel("Boxplots of Individual Signatures",
                                  shiny::sidebarPanel(
                                    shinyWidgets::pickerInput(
                                      'singbox',
                                      'Signature(s)',
                                      choices = siglist,
                                      options = list('actions-box' = TRUE),
                                      multiple = TRUE, selected = NULL),
                      shiny::hr(),
                      shiny::selectInput(
                        'singboxcovar', 'Covariate', choices = NULL
                        ),
                      shiny::actionButton('singboxplot', "Plot Boxplot(s)")
                    ),
                    shiny::mainPanel(
                      shiny::plotOutput("boxplotind", height = 500)
                      )
                    ),
                  shiny::tabPanel(
                    "Compare Scoring Methods for a Single Signature",
                    shiny::sidebarPanel(
                      shiny::selectInput('singcomp',
                                         'Select Signature',
                                         choices = siglist),
                      shiny::selectInput('compassay',
                                         'Select Assay',
                                         choices = NULL),
                      shiny::selectInput('compcovar',
                                         "Select Covariate",
                                         choices = NULL),
                      shinyWidgets::pickerInput('compalg',
                                                'Select Algorithms',
                                                choices = c(
                                                  "GSVA", "ssGSEA",
                                                  "singscore", "PLAGE",
                                                  "Zscore", "ASSIGN"),
                                                options = list(
                                                  'actions-box' = TRUE),
                                                multiple = TRUE,
                                                selected = NULL),
                             shiny::actionButton(
                               'compplot',
                               "Plot Heatmap")
                      ),
                    shiny::mainPanel(
                      shiny::plotOutput("heatcomp", height = 400)
                      )
                    )
                  )
