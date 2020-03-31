observe({
  if(is.null(vals$profilerdat)){
    updatePickerInput(session, 'selectsigs', choices=NULL)
  }
  else{isolate(updatePickerInput(session, 'selectsigs', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))), selected=subset(siglist, siglist %in% colnames(colData(vals$profilerdat)))))
  }
  })

observe({
  updateSelectInput(session, 'allheatcovar', choices = vals$covars)
})

observeEvent(input$allheatplot,
             {output$allheat=renderPlot({
               colors <- RColorBrewer::brewer.pal(6, "Spectral")
               col.me <- circlize::colorRamp2(seq(from = -2, to = 2,
                                                  length.out = 6), colors)

               isolate({signatureHeatmap(vals$profilerdat, name = "Heatmap of Signatures",
                                signatureColNames = names(TBsignatures[input$selectsigs]),
                                annotationColNames = input$allheatcovar,
                                scale = TRUE,
                                showColumnNames = TRUE,
                                choose_color = col.me)})
             })})

observe({
  updateSelectInput(session, 'singheatcovar', choices = vals$covars)
})

observe({
  if(is.null(vals$profilerdat)){
    updatePickerInput(session, 'singheat', choices=NULL)
  }
  else{updatePickerInput(session, 'singheat', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))))
  }
})

observe({
  updatePickerInput(session, 'genes', choices = TBsignatures[input$singheat], selected = TBsignatures[[input$singheat]])
})

observeEvent(input$singheatplot,
             {output$indheat=renderPlot({
               isolate({print(signatureGeneHeatmap(inputData = vals$profilerdat, useAssay = input$profassay,
                                          input$genes,
                                          signatureColNames = input$singheat,
                                          annotationColNames = input$singheatcovar,
                                          showColumnNames = TRUE))})
             })})

observe({
  if(is.null(vals$profilerdat)){
    updatePickerInput(session, 'singbox', choices=NULL)
  }
  else{updatePickerInput(session, 'singbox', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))))
  }
})

observe({
  updateSelectInput(session, "singboxcovar", choices = vals$covars)
})

observeEvent(input$singboxplot,
             {output$boxplotind=renderPlot({
               isolate({print(signatureBoxplot(vals$profilerdat, signatureColNames = input$singbox,
                                      annotationColName = input$singboxcovar))})
             })})

observe({
  if(is.null(vals$profilerdat)){
    updateSelectInput(session, 'singcomp', choices=NULL)
  }
  else{updateSelectInput(session, 'singcomp', choices = subset(siglist, siglist %in% colnames(colData(vals$profilerdat))))
  }
})

observe({
  updateSelectInput(session, "compassay", choices = vals$datassays)
})

observe({
  updateSelectInput(session, "compcovar", choices = vals$covars)
})

observeEvent(input$compplot,
             {output$heatcomp=renderPlot({
               isolate({suppressWarnings(compareAlgs(vals$profilerdat, annotationColName = input$compcovar,
                                            scale = TRUE,
                                            algorithm = input$compalg,
                                            useAssay = input$compassay,
                                            signatures = TBsignatures[input$singcomp],
                                            show.pb = TRUE))})
             })})
