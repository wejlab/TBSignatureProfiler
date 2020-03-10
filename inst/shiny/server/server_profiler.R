output$covars=renderTable({
  vals$covars
}, colnames = F)

output$assays=renderTable({
  vals$datassays
}, colnames = F)

observeEvent(input$runprofiler,{
  output$visdat=DT::renderDataTable(
    head(as.data.frame(colData(vals$profilerdat))), options=list(scrollX=T)
  )
  })

observeEvent(input$mkassay,{
  if(input$newassay=='log'){
    vals$tbdat=mkAssay(vals$tbdat, input_name = 'counts', log = T, counts_to_CPM = F)
    vals$datassays=names(assays(vals$tbdat))
  }
  else if(input$newassay=='cpm'){
    vals$tbdat=mkAssay(vals$tbdat, input_name = 'counts')
    vals$datassays=names(assays(vals$tbdat))
  }
  else if(input$newassay=='logcpm'){
    vals$tbdat=mkAssay(vals$tbdat, input_name = 'counts', log = T)
    vals$datassays=names(assays(vals$tbdat))
  }
})

observe({
  updateSelectInput(session, 'profassay', choices = vals$datassays)
})

observeEvent(input$runprofiler,{
  vals$profilerdat=isolate(runTBsigProfiler(vals$tbdat, useAssay = input$profassay, signatures = TBsignatures, algorithm = input$profalg, combineSigAndAlgorithm = T, parallel.sz = 4))
    })


