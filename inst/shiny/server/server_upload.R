shiny::observeEvent(input$upload, {
  if (input$dat == 'samp'){
    if (input$sampdat == 'tbhiv'){
      vals$tbdat <- TB_hiv
      vals$coldat <- colData(vals$tbdat)
      vals$covars <- colnames(colData(vals$tbdat))
      vals$datassays <- names(assays(vals$tbdat))
    } else if (input$sampdat == 'tbind'){
      vals$tbdat <- TB_indian
      vals$coldat <- colData(TB_indian)
      vals$covars <- colnames(colData(vals$tbdat))
      vals$datassays <- names(assays(vals$tbdat))
      }
    }
  else if (input$dat == 'updat'){
    if (input$updattype == 'rawdat'){
      countdat <- read.table(input$countsfile$datapath, header = TRUE,
                            row.names = 1)
      metadat <- read.table(input$metdatfile$datapath, header = TRUE)
      vals$tbdat <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(countdat)),
                                        colData = metadat)
      vals$coldat <- SummarizedExperiment::colData(vals$tbdat)
      vals$covars <- colnames(SummarizedExperiment::colData(vals$tbdat))
      vals$datassays <- names(SummarizedExperiment::assays(vals$tbdat))
      } else if (input$updattype == 'sumexp'){
        tb1 <- shiny::get(load(input$rdsFile$datapath))
        vals$tbdat <- print(tb1)
        vals$coldat <- SummarizedExperiment::colData(vals$tbdat)
        vals$covars <- colnames(SummarizedExperiment::colData(vals$tbdat))
        vals$datassays <- names(SummarizedExperiment::assays(vals$tbdat))
      }
    }
  })

shiny::observeEvent(input$upload, {
  output$counttable <- DT::renderDataTable(
    as.data.frame(SummarizedExperiment::assays(head(vals$tbdat, n = 20L))),
    options = list(scrollX = TRUE)
    )
  vals$profilerdat <- NULL
  })

shiny::observeEvent(input$upload, {
  output$visdat <- renderDataTable(NULL)
  output$allheat <- NULL
  output$indheat <- NULL
  output$boxplotind <- NULL
  output$heatcomp <- NULL
  output$bootbox <- NULL
  output$rocsep <- NULL
})
