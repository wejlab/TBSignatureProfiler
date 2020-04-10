shiny::observe({
  if (is.null(vals$profilerdat)) {
    shinyWidgets::updatePickerInput(session, 'selectsigs', choices = NULL)
    } else{shiny::isolate(shinyWidgets::updatePickerInput(session,
                                                          'selectsigs',
                                                          choices = subset(
                                                            siglist,
                                                            siglist %in% colnames(
                                                              SummarizedExperiment::colData(
                                                                vals$profilerdat))),
                                                          selected = subset(siglist,
                                                                            siglist %in% colnames(
                                                                              SummarizedExperiment::colData(
                                                                                vals$profilerdat)))))
      }
  })

shiny::observe({
  shiny::updateSelectInput(session, 'allheatcovar', choices = vals$covars)
})

shiny::observeEvent(input$allheatplot, {
  output$allheat <- shiny::renderPlot({
    colors <- RColorBrewer::brewer.pal(6, "Spectral")
    col.me <- circlize::colorRamp2(base::seq(from = -2, to = 2,
                                             length.out = 6), colors)

    shiny::isolate({signatureHeatmap(vals$profilerdat,
                                     name = "Heatmap of Signatures",
                                     signatureColNames = names(
                                       TBsignatures[input$selectsigs]),
                                     annotationColNames = input$allheatcovar,
                                     scale = TRUE,
                                     showColumnNames = TRUE,
                                     choose_color = col.me)
      })
    })
  })

shiny::observe({
  shiny::updateSelectInput(session, 'singheatcovar', choices = vals$covars)
})

shiny::observe({
  if (is.null(vals$profilerdat)){
    shinyWidgets::updatePickerInput(session, 'singheat', choices = NULL)
  } else{shinyWidgets::updatePickerInput(session, 'singheat',
                                         choices = subset(siglist,
                                                          siglist %in% colnames(
                                                            SummarizedExperiment::colData(
                                                              vals$profilerdat))))
    }
  })

shiny::observe({
  shinyWidgets::updatePickerInput(session, 'genes',
                                  choices = TBsignatures[input$singheat],
                                  selected = TBsignatures[[input$singheat]])
  })

shiny::observeEvent(input$singheatplot, {
  output$indheat <- shiny::renderPlot({
    shiny::isolate({
      print(signatureGeneHeatmap(inputData = vals$profilerdat,
                                 useAssay = input$profassay,
                                 input$genes,
                                 signatureColNames = input$singheat,
                                 annotationColNames = input$singheatcovar,
                                 showColumnNames = TRUE))
      })
    })
  })

shiny::observe({
  if (is.null(vals$profilerdat)){
    shinyWidgets::updatePickerInput(session, 'singbox', choices = NULL)
    } else {
      shinyWidgets::updatePickerInput(session, 'singbox',
                                      choices = subset(siglist,
                                                       siglist %in% colnames(
                                                         colData(vals$profilerdat))))
      }
  })

shiny::observe({
  shinyWidgets::updateSelectInput(session, "singboxcovar",
                                  choices = vals$covars)
})

shiny::observeEvent(input$singboxplot, {
  output$boxplotind <- shiny::renderPlot({
    shiny::isolate({print(signatureBoxplot(vals$profilerdat,
                                           signatureColNames = input$singbox,
                                           annotationColName = input$singboxcovar))})
    })
  })

shiny::observe({
  if (is.null(vals$profilerdat)){
    shiny::updateSelectInput(session, 'singcomp', choices = NULL)
    } else{shiny::updateSelectInput(session, 'singcomp',
                                    choices = subset(siglist,
                                                     siglist %in% colnames(
                                                       SummarizedExperiment::colData(
                                                         vals$profilerdat))))
      }
  })

shiny::observe({
  shiny::updateSelectInput(session, "compassay", choices = vals$datassays)
})

shiny::observe({
  shiny::updateSelectInput(session, "compcovar", choices = vals$covars)
})

shiny::observeEvent(input$compplot, {
  output$heatcomp <- shiny::renderPlot({
    shiny::isolate({suppressWarnings(compareAlgs(vals$profilerdat,
                                                 annotationColName = input$compcovar,
                                                 scale = TRUE,
                                                 algorithm = input$compalg,
                                                 useAssay = input$compassay,
                                                 signatures = TBsignatures[input$singcomp],
                                                 show.pb = TRUE))
      })
    })
  })
