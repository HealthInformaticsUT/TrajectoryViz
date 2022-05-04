#' The main function running the Shiny App
#'
#' @param inputData THE INPUT DATA as read_csv("my_data.csv")
#'
#' @export
trajectoryViz <- function(inputData) {
  lapply(c("readr", "tidyr", "dplyr", "ggplot2",
           "shinydashboard", "shiny", "shinyjs",
           "sunburstR", "ggiraph", "stringr", "tibble", "DT", "devtools"),
         require, character.only = TRUE)

  dataTables <- trajectoryDataPrep(inputData)
  freqPaths <- dataTables[[1]]
  patPaths <- dataTables[[2]]
  patStateLevel <- dataTables[[3]]
  colorsDef <- dataTables[[4]]
  labels <- dataTables[[5]]
  colors <- dataTables[[6]]
  ###SHINY UI### -----
  # Header ------------
  header <- dashboardHeader()
  # Sidebar ------------
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Plots",
               tabName = "Plots",
               icon = icon("chart-line")
      ),
      menuItem("Tables",
               tabName = "tables",
               icon = icon("table"),
               badgeLabel = "Data!",
               badgeColor = "light-blue"
      ),
      menuItem("Info",
               tabName = "info",
               icon = icon("question")
      )
    )
  )
  # Body ------------
  body <- dashboardBody(
    fluidRow(
      useShinyjs(),
      tabItems(
        tabItem( # tab item 1 start
          tabName = "Plots",
          fluidRow(
            box(
              width = 4,
              title = "Sequences on a Sunburst Chart",
              sunburstOutput("sunburst1", width = "100%")

            ),
            box(
              width = 4,
              title = "Out of Cohorts excluded",
              girafeOutput("clustPlots")
            ),
            box(
              width = 4,
              title = sprintf("Aligned by selected STATE start"),
              girafeOutput("drugLevel")
            ),

            # ROW
            box(
              width = 4,
              title = "All patients on that path",
              plotOutput("pats1", brush = brushOpts(
                id = "pats1_brush",
                resetOnNew = TRUE
              ))
            ),
            box(
              width = 4,
              title = "Out of Cohorts included",
              girafeOutput("clustPlotsO")
            ),
            box(
              width = 4,
              title = "Aligned by selected STATE start",
              girafeOutput("drugLevel0")
            ),

            ### Row - numbers
            valueBox(
              width = 4,
              value = textOutput("stats1"),
              subtitle = "STATE 1: Average(median) length",
              color = "navy"
            ),
            valueBox(
              width = 4,
              value = textOutput("stats2"),
              subtitle = "STATE 2: Average(median) length",
              color = "navy"
            ),
            valueBox(
              width = 4,
              value = textOutput("stats3"),
              subtitle = "STATE 3: Average(median) length",
              color = "navy"
            )
          ) # Fluidrow end
        ), # TabItem end
        tabItem( # tab item 1 start
          tabName = "tables",
          # fluidRow(
          tabBox(
            width = 12,
            id = "tabset1",
            tabPanel("Frequences",
                     downloadButton(
                       outputId = "freq_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tableFreq")),
            tabPanel("Paths",
                     downloadButton(
                       outputId = "paths_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tablePaths")),
            tabPanel("Levels",
                     downloadButton(
                       outputId = "levels_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tableLevels")),
            tabPanel("Source",
                     downloadButton(
                       outputId = "source_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tableAll"))
          ) # tabBox end
          # )# Fluidrow end
        ), # TabItem end
        tabItem( # tab item 1 start
          tabName = "info",
          tabsetPanel(type = "tabs",
                      tabPanel(title = sprintf("About Plots tab"),
                               imageOutput("plotImage"),
                               textOutput("infoPlots")),
                      tabPanel(title = sprintf("About Tables tab"),
                               tags$style(
                                 "p, div { color: #444444;
                                  padding-left: 1px;
                                 }"
                               ),
                               h3("Frequences"),
                               div("A block"),
                               h3("Paths"),
                               div("A block"),
                               h3("Levels"),
                               div("A block"),
                               h3("Source"),
                               div("A block"),
                               textOutput("infoTables"))
          )
        ) # TabItem end
      ) # tabItems end
    ) # fluidRow end
  ) #  dashboardBody end

  # UI ----
  ui <- dashboardPage(
    skin = "purple",
    header,
    sidebar,
    body)

  # Server ----

  server <- function(input, output) {
    #data----

    pathList <- reactive({input$sunburst1_click})
    aSelected <- reactive({input$clustPlots0_selected})
    bSelected <- reactive({input$clustPlots_selected})
    output$aSelected <- renderText(aSelected())
    output$bSelected <- renderText(bSelected())
    output$path <- renderText(pathList())

    output$sunburst1 <- renderSunburst({
      add_shiny(sunburst(freqPaths,
                         count = TRUE,
                         colors = list(range = c(colors, "#cccccc","#cccccc"), domain = c(labels, "OUT OF COHORT", "End")),
                         legend = list(w=  100, h = 20, s = 5), #TRUE,
                         breadcrumb = htmlwidgets::JS(("function(x) {return x;}"))
      ))
    })

    output$clustPlotsO <- renderGirafe({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        clustPlotsO(patStateLevel, pathSelected(), patPaths, colorsDef)
      })
    })

    output$clustPlots <- renderGirafe({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        clustPlots(patStateLevel, pathSelected(), patPaths, colorsDef)
      })
    })

    output$drugLevel0 <- renderGirafe ({
      if (is.null(input$clustPlots_selected)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        drugLevel0(patStateLevel, pathSelected(), patPaths, colorsDef, bSelected())
      })
    }) ## end drugLevel0

    output$drugLevel <- renderGirafe ({
      if (is.null(input$clustPlots_selected)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        drugLevel(patStateLevel, pathSelected(), patPaths, colorsDef, bSelected())
      })
    }) ## end drugLevel

    # observeEvent(input$clustPlots0_selected, {
    output$plotSort0 <- renderGirafe ({
      if (is.null(input$clustPlots_selected)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        plotSort0(patStateLevel, pathSelected(), patPaths, colorsDef, bSelected())
      })
    })#}
    # ) # end observeEvent clustPlots0_selected

    #observeEvent(input$clustPlots_selected, {
    output$plotSort <- renderGirafe ({
      if (is.null(input$clustPlots_selected)) {
        return(NULL)
      }
      isolate({
        pathSelected <- renderText(pathList())
        plotSort(patStateLevel, pathSelected(), patPaths, colorsDef, bSelected())
      })
    })#}
    # ) # end observeEvent clustPlots_selected

    ## PATHS by total LENGTH ----
    output$pats1 <- renderPlot({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        selectPatientsStartAlpha(patStateLevel, pathSelected(), patPaths, colorsDef)
      })
    })
    ## ZOOM in CLUSTERS ----
    ranges1 <- reactiveValues(x = NULL, y = NULL)
    output$clustZoom <- renderPlot({
      if (is.null(input$clustZoom_brush)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        ranges <- ranges1
        clustPlots(patStateLevel, pathSelected(), patPaths, colorsDef) +
          coord_cartesian(xlim = ranges1$x, ylim = ranges1$y, expand = FALSE)
      })
    })
    observe({
      brush <- input$clustZoom_brush
      if (!is.null(brush)) {
        ranges1$x <- c(brush$xmin, brush$xmax)
        ranges1$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges1$x <- NULL
        ranges1$y <- NULL
      }
    })
    ## ZOOM IN PATIENTS ALL ----
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    output$pats1zoom <- renderPlot({
      if (is.null(input$pats1_brush)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        ranges <- ranges2
        selectPatientsStartAlpha(patStateLevel, pathSelected(), patPaths, colorsDef) +
          coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = FALSE)
      })
    })
    observe({
      brush <- input$pats1_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })

    ## Shiny Infoblocks----
    output$stats1 <- renderText({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        stats1(pathSelected(), patPaths)
      })
    })
    output$stats2 <- renderText({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        stats2(pathSelected(), patPaths)
      })
    })
    output$stats3 <- renderText({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        stats3(pathSelected(), patPaths)
      })
    })
    output$stats4 <- renderText({
      if (is.null(input$sunburst1_click)) {
        return(NULL)
      }
      isolate({
        pathList <- reactive({
          input$sunburst1_click
        })
        pathSelected <- renderText(pathList())
        stats4(pathSelected(), patPaths)
      })
    })
    ## Shiny Tables Downloads ----
    tablePaths <- reactive({patPaths})
    tableAll <- reactive({inputData})
    tableFreq <- reactive({freqPaths})
    tableLevels <- reactive({patStateLevel})

    ## Shiny Table Outputs ----
    output$tableAll <- DT::renderDataTable({
      datatable(inputData, filter = "top", options = list(scrollX = TRUE))
    })
    output$tablePaths <- DT::renderDataTable({
      datatable(patPaths, width ="95%",  height = "auto", filter = "top", options = list(scrollX = TRUE))
    })
    output$tableLevels <- DT::renderDataTable({
      datatable(patStateLevel, width ="95%",  height = "auto", filter = "top", options = list(scrollX = TRUE))
    })
    output$tableFreq <- DT::renderDataTable({
      tableFreq <- freqPaths
      datatable(tableFreq, filter = "top", options = list(scrollX = TRUE)) %>%
        formatStyle(
          "freq",
          background = styleColorBar(c(0, max(tableFreq$freq)), "lightblue")
        )
    })

    output$paths_download_button <- shiny::downloadHandler(
      filename = paste0("paths-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tablePaths(), file)
      }
    )
    output$source_download_button <- shiny::downloadHandler(
      filename = paste0("source-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tableAll(), file)
      }
    )
    output$freq_download_button <- shiny::downloadHandler(
      filename = paste0("freq-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tableFreq(), file)
      }
    )
    output$levels_download_button <- shiny::downloadHandler(
      filename = paste0("levels-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tableLevels(), file)
      }
    )
    ## Shiny INFO tab ----
    output$infoPlots <- renderText("")
    output$infoTables <- renderText("")
    output$plotImage <- renderImage(list(src='plot-ov1.png', width ="95%",  height = "auto"),
                                    deleteFile = FALSE)
  }

  shiny::shinyApp(ui, server)
}

