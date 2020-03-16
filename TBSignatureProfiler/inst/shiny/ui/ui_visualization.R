navbarMenu("Signature Visualization",
           tabPanel("Heatmap with Selected Signatures",
                    sidebarPanel(
                      pickerInput('selectsigs', label='Select Signatures for Profiler',choices=siglist, options=list('actions-box'=T),multiple=T, selected = siglist),
                      selectInput('allheatcovar', "Select Covariate", choices = NULL),
                      actionButton('allheatplot', "Plot Heatmap")
                    ),
                    mainPanel(
                      plotOutput("allheat", height = 800)
                    )
           ),
           tabPanel("Heatmaps of Individual Signatures",
                    sidebarPanel(
                      pickerInput('singheat', 'Signature', choices =siglist,multiple=F, selected =NULL),
                      pickerInput('genes', label='Genes',choices=NULL, options=list('actions-box'=T),multiple=T),
                      hr(),
                      selectInput('singheatcovar', 'Covariate', choices =NULL),
                      actionButton('singheatplot', "Plot Heatmap(s)")
                    ),
                    mainPanel(
                      plotOutput("indheat", height = 800)
                    )
           ),
           tabPanel("Boxplots of Individual Signatures",
                    sidebarPanel(
                      pickerInput('singbox', 'Signature(s)', choices =siglist, options=list('actions-box'=T),multiple=T, selected =NULL),
                      hr(),
                      selectInput('singboxcovar', 'Covariate', choices =NULL),
                      actionButton('singboxplot', "Plot Boxplot(s)")
                    ),

                    mainPanel(
                      plotOutput("boxplotind", height = 500)
                    )
           ),
           tabPanel("Compare Scoring Methods for a Single Signature",
                    sidebarPanel(
                      selectInput('singcomp', 'Select Signature', choices =siglist ),
                      selectInput('compassay', 'Select Assay', choices = NULL),
                      selectInput('compcovar', "Select Covariate", choices = NULL),
                      pickerInput('compalg', 'Select Algorithms', choices = c("GSVA", "ssGSEA", "singscore", "PLAGE", "Zscore", "ASSIGN"),options=list('actions-box'=T),multiple=T, selected =NULL),
                      actionButton('compplot', "Plot Heatmap")
                    ),
                    mainPanel(
                      plotOutput("heatcomp", height = 400)
                    )
           )
)
