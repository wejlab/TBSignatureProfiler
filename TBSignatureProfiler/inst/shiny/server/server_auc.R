output$ttest=DT::renderDataTable({
  tableAUC(vals$profilerdat,
           annotationColName = "Disease",
           signatureColNames = names(TBsignatures),
           num.boot = 100,
           pb.show = FALSE)
})

observe({
  if(is.null(vals$profilerdat)){
    updatePickerInput(session, 'bootsigs', choices=NULL)
  }
  else{updatePickerInput(session, 'bootsigs', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))))
  }
})

observeEvent(input$bootplot,{
  colors.box=RColorBrewer::brewer.pal(12,"Set3")
  colors.rep=rep(colors.box,3)
  l.sig=length(input$bootsigs)
  output$bootbox=renderPlot({
    isolate({compareBoxplots(vals$profilerdat, annotationColName = "Disease",
                    signatureColNames = names(TBsignatures[input$bootsigs]),
                    pb.show = FALSE, fill.col = colors.rep[1:l.sig], rotateLabels = T)})
  })
})

observe({
  if(is.null(vals$profilerdat)){
    updatePickerInput(session, 'singroc', choices=NULL)
  }
  else{updatePickerInput(session, 'singroc', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))))
  }
})

observeEvent(input$rocplot,{
  output$rocsep=renderPlot({
    isolate({print(signatureROCplot_CI(inputData = vals$profilerdat,
                              signatureColNames = input$singroc,
                              annotationColName = "Disease",
                              name = paste("ROC plot")))})
  })
})
