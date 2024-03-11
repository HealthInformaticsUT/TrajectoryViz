#' The main function running the Shiny App
#'
#' @param inputData THE INPUT DATA as read_csv("dataASCH.csv")
#'
#' @export
#'
trajectoryViz <- function() { ###
  lapply(c("readr", "tidyr", "tidyverse", "dplyr", "ggplot2",
           "shinydashboard", "shinydashboardPlus", "shiny", "shinyjs",
           "sunburstR", "ggiraph", "stringr", "tibble", "DT", "devtools",
           "ggrepel", "shinycssloaders", "plotly"),
         require, character.only = TRUE)
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
    ),
    hr(),
    br(),
    collapsed = TRUE
  )
  # Body ------------
  body <- dashboardBody(
    fluidRow(
      useShinyjs(),
      tabItems(
        tabItem( ##### tab item 1 start #####
          tabName = "Plots",
          fluidRow(
            hr(),
            #### ROW1 ####
             box(
               solidHeader = FALSE,
               title = "Data overview",
               background = NULL,
               width = 6,
               status = "info",
               footer = fluidRow(
                 column(
                   width = 6,
                   fileInput(inputId = "user_uploaded",
                             label = "Choose CSV File for input data",
                             accept = ".csv", width = 400,
                             buttonLabel = "Browse ...",
                             placeholder = "No file selected",
                             multiple = FALSE),
                   valueBox(
                     width = 12,
                     value = textOutput("n_patients1"),
                     subtitle = "Patients",
                     color = "light-blue"
                   )
                 ),
                 column(
                   width = 6,
                   selectInput("startingState",
                               label= "Select starting state:",
                               choices = c("Start as imported"),
                               selected = "Start as imported"),
                   br(),
                   valueBox(
                     width = 12,
                     value = textOutput("n_traj1"),
                     subtitle = "Unique sequences",
                     color = "light-blue"
                   )
                 )
               )
             ),
            box(
              width = 6,
              title = "Sequences Divided on a Sunburst Chart",
              withSpinner(sunburstOutput("sunburst2", width = "100%"), color = "#3D8FBE", type = 7, size=0.6)
            ),
            box(
              width = 12
            ),
            #### ROW2 ####
            box(
              width = 6,
              title = "Patient-level Sequences",
              withSpinner(girafeOutput("clustPlots2"), color = "#3D8FBE", type = 7, size=0.6)
            ),
            box(
              width = 6,
              title = "Patient-level Sequences Aligned and Sorted",
              withSpinner(girafeOutput("drugLevel02"), color = "#3D8FBE", type = 7, size=0.6),
              column(
                width = 4,
                selectInput("arrBy2",
                            choices = c(NULL),
                            selected = NULL,
                            p("Arrange by distance from:"),
                            width = "200px")
              ),
              column(
                width = 4,
                numericInput("vline1", p("Indicate distance"), min= 0, 365, width = "200px")
                 ),
              column(
                width = 4,
                radioButtons("befAft", p("From alignemnt point"), c("Before", "After", "Before or After"), width = "200px")
              )
            ),

            box(
              width = 12
            ),
            box(
              width = 6,
              title = "Funnel",
              plotlyOutput("funnel")
            ),
            box(
              width = 12
            ),
          ) # Fluidrow end
        ), # TabItem end
        tabItem( # tab item 2 start #####
          tabName = "tables",
          # fluidRow(
          tabBox(
            width = 12,
            id = "tabset1",
            tabPanel("Unique Path Frequences",
                     downloadButton(
                       outputId = "freq2_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tableFreqPercent")),
            tabPanel("Individual Paths",
                     downloadButton(
                       outputId = "paths_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tablePaths")),
            tabPanel("Levels of All Paths",
                     downloadButton(
                       outputId = "levels_download_button",
                       label = "Download CSV"
                     ),
                     DT::dataTableOutput("tableLevels")),
            tabPanel("Source Data",
                     DT::dataTableOutput("tableAll1")),
          ) # tabBox end
          # )# Fluidrow end
        ), # TabItem end
        tabItem( # tab item 3 start #####
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
                               h3("Unique Path Frequences"),
                               div("Input for the sunburst. Each row is a unique sequence with its frequency and % in the dataset"),
                               h3("Individual Paths"),
                               div("For the patient-level trajectory plot. Each row indicates a patient and its treatment trajectory details."),
                               div("In columns: ID, sequence as a string, each STATE (L1 - L11), length of each state (Len1 - Len11)"),
                               h3("Levels of All Paths"),
                               div("For the patient-level trajectory plot."),
                               div("Each row indicates a discrete period (length pre-defined in Cohort2Trajectory) from the patient's treatment trajectory."),
                               h3("Source Data"),
                               div("Input data that was used to run trajectoryViz()."),
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

  server <- function(input, output, session) {

    ## Read in uploaded data file ####
    #### Default limit is 5MB a file. If needed to resize, run in console before executing the Shiny:
    #### options(shiny.maxRequestSize = 20 * 1024^2)
    inputData1 <- reactive({
      req(input$user_uploaded)
      inputData <- read_csv(input$user_uploaded$datapath, col_types = cols())
      if("STATE_LABEL" %in% colnames(inputData)){
        inputData = inputData %>%
          rename("STATE" = "STATE_LABEL")
      }
      return(inputData)
    })

    ## Data preparation from Uploaded file ####
    dataTables1 <- reactive({
      if (input$startingState == "Start as imported")({
        return(trajectoryDataPrep(input = inputData1()))
        })
      isolate(
      if (input$startingState != "Start as imported")({
          startState <- reactive({input$startingState})
          x <- reactive(trajectoryDataPrep(input = inputData1()))
          patStateLevel <- reactive({x()[[3]]})
          inputDataCut <- reactive(cutStart(patStateLevel(), startState()))
          return(trajectoryDataPrep(input = inputDataCut()))})
      )
      })

    freqPaths1 <- reactive(dataTables1()[[1]])
    patPaths1 <- reactive(dataTables1()[[2]])
    patStateLevel1 <- reactive({dataTables1()[[3]]})
    freqPathsPercent1 <- reactive(dataTables1()[[7]])


    labels <- reactive(c(unique(inputData1()$STATE), "End")) # "End" is required for RSunburst package input
    labels1 <- reactive(labels()[!labels() %in% c("START", "EXIT", "OUT OF COHORT", "End")])
    colors_all <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73",
                    "#0072B2", "#D55E00", "#CC79A7", "#004949",
                    "#b66dff", "#924900", "#ff6db6", "#490092", "#601A4A" )
    colors1 <- reactive(colors_all[1:length(labels1())])
    colorsDef1 <- reactive(setNames(c(colors1(),"#cccccc"), c(labels1(), "OUT OF COHORT")))






    n_patients1 <- reactive({
      data <- patStateLevel1()
      n_distinct(data$SUBJECT_ID)
    })
    n_traj1 <- reactive({
      data <- freqPaths1()
      nrow(data)
    })
    output$n_patients1 <- renderText({return(n_patients1())})
    output$n_traj1 <- renderText({return(n_traj1())})

    ## Update STATES for dropDown selections ####
    all_states <- reactive({
      patStateLevel1 <- patStateLevel1()
      unique(str_c(patStateLevel1$STATE, "-", patStateLevel1$SEQ_ORDINAL))
    })
    arrByList <- reactive({
      patStateLevel1 <- patStateLevel1()
      unique(str_c(patStateLevel1$STATE, "-", patStateLevel1$SEQ_ORDINAL))
    })
    observeEvent(input$user_uploaded, updateSelectInput(session, "startingState", choices = c("Start as imported", all_states())))
    observeEvent(input$user_uploaded, updateSelectInput(session, "arrBy2", choices = c(arrByList())))

    ## Reactive user input data ####
    startState <- reactive({input$startingState})
    arrBy2 <- reactive({input$arrBy2})
    vline1 <- reactive({input$vline1})
    befAft <- reactive({input$befAft})
    pathList2 <- reactive({input$sunburst2_click})
    bSelected2 <- reactive({input$clustPlots2_selected})
    output$bSelected2 <- renderText(bSelected2())
    output$path2 <- renderText(pathList2())

    ## PLOT Outputs #####
    output$sunburst2 <- renderSunburst({
      validate(
        need((!is.null(input$user_uploaded)), "Upload a data set.")
      )
          return(add_shiny(
            htmlwidgets::onRender(
              sunburst(freqPaths1(),
                       count = TRUE,
                       colors = list(range = c(colors1(), "#cccccc","#cccccc"), domain = c(labels1(), "OUT OF COHORT", "End")),
                       legend = list(w=  200, h = 20, s = 5), #TRUE,
                       breadcrumb = htmlwidgets::JS(("function(x) {return x;}"))),
              "
    function(el,x){
    d3.select(el).select('.sunburst-togglelegend').property('checked', true);
    d3.select(el).select('.sunburst-legend').style('visibility', '');
    }
    "
            )
          )
          )
    })

    output$clustPlots2 <- renderGirafe({
      validate(
        need((!is.null(input$user_uploaded)), "Upload a data set.")
      )
          pathSelected <- renderText(str_c(pathList2(), collapse = "---"))
          return(clustPlots(patStateLevel1(), pathSelected(), patPaths1(), colorsDef1()))
    })

    output$drugLevel02 <- renderGirafe ({
      pathSelected <- renderText(str_c(pathList2(), collapse = "---"))
        validate(
          need((!is.null(input$user_uploaded)), "Upload a data set."),
          need((!is.null(input$clustPlots2_selected)), "Choose the State to align by on the plot left."),
          need((!is.null(patPaths1())), "No data to be displayed. Filter again.")
        )
        return(drugLevel0(patStateLevel1(), pathSelected(), patPaths1(), colorsDef1(), bSelected2(), arrBy2(), vline1(), befAft()))
    })

    output$funnel <- renderPlotly ({
      validate(
        need((!is.null(input$user_uploaded)), "Upload a data set."),
        need((!is.null(input$clustPlots2_selected)), "Choose the State to align by on the plot left.")
      )
      pathSelected <- renderText(str_c(pathList2(), collapse = "---"))
      return(funnel(patStateLevel1(), pathSelected(), patPaths1(), bSelected2(), arrBy2(), vline1(), befAft()))
    })

    ## Table Outputs ----
    output$tableFreqPercent <- DT::renderDataTable({
      validate(
        need((!is.null(input$user_uploaded)), "No data set uploaded yet."),
      )
      tableFreqPercent <- freqPathsPercent1()
      datatable(freqPathsPercent1(), filter = "top", options = list(scrollX = TRUE)) %>%
        formatStyle("freq", background = styleColorBar(c(0, max(tableFreqPercent$freq)), "lightblue"))
    })
    output$tableAll1 <- DT::renderDataTable({
      validate(
        need((!is.null(input$user_uploaded)), "No data set uploaded yet."),
      )
      datatable(inputData1(), filter = "top", options = list(scrollX = TRUE))
    })
    output$tablePaths <- DT::renderDataTable({
      validate(
        need((!is.null(input$user_uploaded)), "No data set uploaded yet."),
      )
      datatable(patPaths1(), width ="95%",  height = "auto", filter = "top", options = list(scrollX = TRUE))
    })
    output$tableLevels <- DT::renderDataTable({
      validate(
        need((!is.null(input$user_uploaded)), "No data set uploaded yet."),
      )
      datatable(patStateLevel1(), width ="95%",  height = "auto", filter = "top", options = list(scrollX = TRUE))
    })

    ## Tables Downloads ----
    tableFreqPercent <- reactive({freqPathsPercent1})
    tableAll <- reactive({inputData1})
    tablePaths <- reactive({patPaths1})
    tableLevels <- reactive({patStateLevel1})

    ## Table Download Buttons ####
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
    output$freq2_download_button <- shiny::downloadHandler(
      filename = paste0("freq2-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tableFreqPercent(), file)
      }
    )
    output$levels_download_button <- shiny::downloadHandler(
      filename = paste0("levels-", Sys.Date(), ".csv"),
      content = function(file) {
        write.csv(tableLevels(), file)
      }
    )
    ## INFO tab ----
    output$infoPlots <- renderText("")
    output$infoTables <- renderText("")
    output$plotImage <- renderImage(list(src='plot-1.png', width ="95%",  height = "auto"),
                                    deleteFile = FALSE)
  }
  shiny::shinyApp(ui, server)
}

