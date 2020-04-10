output$covars <- shiny::renderTable({vals$covars}, colnames = FALSE)

output$assays <- shiny::renderTable({vals$datassays}, colnames = FALSE)

shiny::observeEvent(input$runprofiler, {
  output$visdat <- DT::renderDataTable(
    head(as.data.frame(SummarizedExperiment::colData(vals$profilerdat))),
    options = list(scrollX = TRUE)
    )
  })

shiny::observeEvent(input$mkassay, {
  if (input$newassay == 'log'){
    vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts',
                         log = TRUE, counts_to_CPM = FALSE)
    vals$datassays <- names(SummarizedExperiment::assays(vals$tbdat))
    } else if (input$newassay == 'cpm'){
      vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts')
      vals$datassays <- names(SummarizedExperiment::assays(vals$tbdat))
      } else if (input$newassay == 'logcpm'){
        vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts', log = TRUE)
        vals$datassays <- names(SummarizedExperiment::assays(vals$tbdat))
        }
  })

shiny::observe({
  updateSelectInput(session, 'profassay', choices = vals$datassays)
  })

shiny::observeEvent(input$runprofiler, {
  vals$profilerdat <- shiny::isolate(
    runTBsigProfiler(vals$tbdat,
                     useAssay = input$profassay,
                     signatures = TBsignatures,
                     algorithm = input$profalg,
                     combineSigAndAlgorithm = TRUE,
                     parallel.sz = 4))
  })

shiny::observeEvent(input$runprofiler, {
  output$allheat <- NULL
  output$indheat <- NULL
  output$boxplotind <- NULL
  output$heatcomp <- NULL
  output$bootbox <- NULL
  output$rocsep <- NULL
  })
