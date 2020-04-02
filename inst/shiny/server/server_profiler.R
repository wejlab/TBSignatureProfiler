output$covars <- renderTable({
  vals$covars},
  colnames = FALSE)

output$assays <- renderTable({
  vals$datassays},
  colnames = FALSE)

observeEvent(input$runprofiler, {
  output$visdat <- DT::renderDataTable(
    head(as.data.frame(colData(vals$profilerdat))),
    options = list(scrollX = TRUE)
  )
  })

observeEvent(input$mkassay, {
  if (input$newassay == 'log'){
    vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts',
                         log = TRUE, counts_to_CPM = FALSE)
    vals$datassays <- names(assays(vals$tbdat))
  } else if (input$newassay == 'cpm'){
    vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts')
    vals$datassays <- names(assays(vals$tbdat))
  } else if (input$newassay == 'logcpm'){
    vals$tbdat <- mkAssay(vals$tbdat, input_name = 'counts', log = TRUE)
    vals$datassays <- names(assays(vals$tbdat))
  }
})

observe({
  updateSelectInput(session, 'profassay', choices = vals$datassays)
})

observeEvent(input$runprofiler, {
  vals$profilerdat <- isolate(runTBsigProfiler(vals$tbdat,
                                              useAssay = input$profassay,
                                              signatures = TBsignatures,
                                              algorithm = input$profalg,
                                              combineSigAndAlgorithm = TRUE,
                                              parallel.sz = 4))
    })

observeEvent(input$runprofiler, {
  output$allheat <- NULL
  output$indheat <- NULL
  output$boxplotind <- NULL
  output$heatcomp <- NULL
  output$bootbox <- NULL
  output$rocsep <- NULL
})
