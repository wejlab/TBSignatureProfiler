output$ttest <- DT::renderDataTable({
  tableAUC(vals$profilerdat,
           annotationColName = "Disease",
           signatureColNames = names(TBsignatures),
           num.boot = 100,
           pb.show = FALSE)
})

shiny::observe({
  if (is.null(vals$profilerdat)) {
    shinyWidgets::updatePickerInput(session, "bootsigs", choices = NULL)
  } else {
    shinyWidgets::updatePickerInput(
      session, "bootsigs",
      choices = subset(siglist,
                       siglist %in% colnames(
                         SummarizedExperiment::colData(vals$profilerdat))))
  }
})

shiny::observeEvent(input$bootplot, {
  colors_box <- RColorBrewer::brewer.pal(12, "Set3")
  colors_rep <- rep(colors_box, 3)
  l_sig <- length(input$bootsigs)
  output$bootbox <- shiny::renderPlot({
    shiny::isolate({
      compareBoxplots(
        vals$profilerdat, annotationColName = "Disease",
        signatureColNames = names(TBsignatures[input$bootsigs]),
        pb.show = FALSE, fill.col = colors_rep[1:l_sig],
        rotateLabels = TRUE)
    })
  })
})

shiny::observe({
  if (is.null(vals$profilerdat)) {
    shinyWidgets::updatePickerInput(session, "singroc", choices = NULL)
  } else {
    shinyWidgets::updatePickerInput(
      session, "singroc", choices = subset(
        siglist, siglist %in% colnames(
          SummarizedExperiment::colData(vals$profilerdat))))
  }
})

shiny::observeEvent(input$rocplot, {
  output$rocsep <- shiny::renderPlot({
    shiny::isolate({
      print(signatureROCplot_CI(inputData = vals$profilerdat,
                                signatureColNames = input$singroc,
                                annotationColName = "Disease",
                                name = paste("ROC plot")))
    })
  })
})
