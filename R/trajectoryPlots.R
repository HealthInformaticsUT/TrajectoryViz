#' Function cutting data from selected start
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param startState reactive input from user input
cutStart <- function(patStateLevel, startState){
  selectedStartData = patStateLevel %>%
    group_by(SUBJECT_ID) %>%
    mutate(CONTAINS_LEVEL = any(str_c(STATE, "-", SEQ_ORDINAL) == startState)) %>%
    ungroup() %>%
    filter(CONTAINS_LEVEL) %>%
    select(-CONTAINS_LEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(STATE_START_DATE >= STATE_START_DATE[str_c(STATE, "-", SEQ_ORDINAL) == startState])

  return(selectedStartData)
}

#' Function plotting patient trajectories without OUT OF COHORT periods
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
clustPlots <- function(patStateLevel, pathClicked, patPaths, colorsDef) { # FUNCTION clustPlotsO sunburst click -> plot with NO OUT OF COHORT periods
  data <- patStateLevel

  selected <- patPaths %>%
    filter(str_detect(path1, str_c("^", pathClicked)))

  if(nrow(selected) > 1){
    mat_data <- selected %>%
      select(SUBJECT_ID, starts_with("Len"))
    mat <- mat_data %>%
      as.data.frame() %>%
      column_to_rownames("SUBJECT_ID") %>%
      select(1:length(str_split(str_replace(pathClicked, "---End", ""), "---")[[1]])) %>%
      as.matrix()



    hc <- hclust(dist(mat, method = "euclidean"), method = "complete")
    ord1 <- tibble(
      Rank1 = 1:nrow(mat),
      SUBJECT_ID = as.double(rownames(mat)[hc$order])
    )
  }

  else(
    ord1 <- tibble(
      Rank1 = 1,
      SUBJECT_ID = selected$SUBJECT_ID
    )
  )


  data1 <- data %>%
    filter(SUBJECT_ID %in% selected$SUBJECT_ID)
  data <- left_join(data1, ord1, by = join_by(SUBJECT_ID))

  dataEvents <- data %>%
    group_by(SUBJECT_ID) %>%
    mutate(REL_STATE_START_DATE = as.numeric(STATE_START_DATE - min(STATE_START_DATE))) %>%
    mutate(REL_STATE_END_DATE = as.numeric(STATE_END_DATE - min(STATE_START_DATE))) %>%
    ungroup()

  dataOOC <- dataEvents %>%
    group_by(SUBJECT_ID, Rank1) %>%
    summarise(
      REL_STATE_START_DATE = min(REL_STATE_START_DATE),
      REL_STATE_END_DATE = max(REL_STATE_END_DATE)
    ) %>%
    ungroup()

  # Set geom parameters
  plot <- ggplot(dataEvents) +
    geom_rect(aes(xmin = REL_STATE_START_DATE, xmax = REL_STATE_END_DATE, ymin = Rank1, ymax = Rank1 + 1), fill = colorsDef["OUT OF COHORT"], data = dataOOC) +
    geom_rect_interactive(aes(xmin = REL_STATE_START_DATE, xmax = REL_STATE_END_DATE, ymin = Rank1, ymax = Rank1 + 1, tooltip = str_c(STATE, "-", SEQ_ORDINAL), data_id = str_c(STATE, "-", SEQ_ORDINAL), fill = STATE), data = dataEvents) +
    theme_classic() +
    scale_fill_manual(values = colorsDef) +
    labs(title = sprintf("Selected: %s ... N = %i ", pathClicked, length(selected$SUBJECT_ID))) +
    theme(
      axis.text.x = element_text(face = "plain", size = 8, colour = "black"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(face = "plain", size = 8, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous("Time (days)") +
    scale_y_continuous("Sequences sorted by eucledian distance based on STATE lengths")+
    guides(fill = guide_legend(
      title = "STATE",
      title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
      label.theme = element_text(face = "plain", size = 8, colour = "darkgrey")
    ))

  girafe(ggobj = plot,
         height_svg = 6,
         width_svg = 8,
         options = list(
           opts_tooltip(
             opacity = .8,
             css = "background-color:yellow; color:black; padding:2px; border-radius:5px;"
           ),
           opts_toolbar(position = "topright"),
           opts_selection(type = "single", only_shiny = FALSE, css = "stroke:#fdfd66;stroke-width:1.4"),
           opts_hover(css = "fill:#444444;stroke:#444444;cursor:pointer;")
         ))
}


# id_new = trajectoryDataPrep2(inputData)
# id2_new = trajectoryDataPrep2(inputData2)
#
# id_old = trajectoryDataPrep(inputData)
#
#
# patStateLevel = id2_new[[3]]
# pathClicked = "PC diagnosis---Low risk"
# patPaths = id2_new[[2]]
# colorsDef = id2_new[[4]]
#
# clustPlots2(patStateLevel, pathClicked, patPaths, colorsDef)


#' Function plotting patient trajectories with OUT OF COHORT periods and aligning the trajectories by selected drug
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug
#' @param arrBy reactive user input
#' @param vline reactive user input

drugLevel0 <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug, arrBy, vline, befAft){ #Alignment plot WITH out of cohorts
  selectedPaths <- patPaths %>%
    filter(str_detect(path1, str_c("^", pathClicked)))

  selectedLevels <- patStateLevel %>%
    filter(SUBJECT_ID %in% selectedPaths$SUBJECT_ID) %>%
    group_by(SUBJECT_ID) %>%
    mutate(REL_STATE_START_DATE = as.numeric(STATE_START_DATE - min(STATE_START_DATE))) %>%
    mutate(REL_STATE_END_DATE = as.numeric(STATE_END_DATE - min(STATE_START_DATE))) %>%
    ungroup()

  alignTime = selectedLevels %>%
    group_by(SUBJECT_ID) %>%
    filter(str_c(STATE, "-", SEQ_ORDINAL) == selectedDrug) %>%
    summarise(pathAlign = min(REL_STATE_START_DATE)) %>%
    ungroup()

  arrangeRank = selectedLevels %>%
    group_by(SUBJECT_ID) %>%
    filter(str_c(STATE, "-", SEQ_ORDINAL) == arrBy) %>%
    summarise(pathStart = min(REL_STATE_START_DATE)) %>%
    ungroup() %>%
    inner_join(alignTime, by = join_by(SUBJECT_ID)) %>%
    mutate(Rank1 = rank(pathAlign - pathStart, ties.method = "first"))


  ### added new - to be checked

  arrangeRank <- arrangeRank %>%
    mutate(distance = pathStart-pathAlign)%>%
    filter(!if_any(everything(), is.na))


  if (befAft == "Before"){
    arrangeRank$vlinePos1 <- -vline
    arrangeRank$vlinePos2 <- 0
  } else if (befAft == "After"){
    arrangeRank$vlinePos1 <- 0
    arrangeRank$vlinePos2 <- vline
  } else {
    arrangeRank$vlinePos1 <- -vline
    arrangeRank$vlinePos2 <- vline
  }
  #arrangeRank$vlinePos <- vline

  validSubjects <- arrangeRank %>%
    filter(distance >= vlinePos1 & distance <= vlinePos2)
  nValid <- length(validSubjects$SUBJECT_ID)
  if (nValid == 0){
    validFirstRank <- 0
    validLastRank <- 0
  } else {
    validFirstRank <- min(validSubjects$Rank1)
    validLastRank <- max(validSubjects$Rank1)
  }
  middleValid <- validFirstRank+(validLastRank - validFirstRank)/2
 # percValid <- nValid/length(unique(dataEvents$SUBJECT_ID))*100
  ### check end


  dataEvents <- inner_join(arrangeRank, selectedLevels, by = join_by(SUBJECT_ID)) %>%
    mutate(REL_STATE_START_DATE = REL_STATE_START_DATE - pathAlign) %>%
    mutate(REL_STATE_END_DATE = REL_STATE_END_DATE - pathAlign)

  dataOOC <- dataEvents %>%
    group_by(SUBJECT_ID, Rank1) %>%
    summarise(
      REL_STATE_START_DATE = min(REL_STATE_START_DATE),
      REL_STATE_END_DATE = max(REL_STATE_END_DATE)
    ) %>%
    ungroup()


  p1 <- ggplot() +
    annotate(geom = "rect", xmin = min(dataEvents$REL_STATE_START_DATE), xmax = max(dataEvents$REL_STATE_END_DATE), ymin = validFirstRank-0.5, ymax = validLastRank,
             fill = "palegreen", alpha = 0.2) +
    geom_rect(aes(xmin = REL_STATE_START_DATE, xmax = REL_STATE_END_DATE, ymin = Rank1 -1, ymax = Rank1), fill = colorsDef["OUT OF COHORT"], data = dataOOC) +
    geom_rect_interactive(aes(xmin = REL_STATE_START_DATE, xmax = REL_STATE_END_DATE, ymin = Rank1 - 1, ymax = Rank1, tooltip = str_c(STATE, "-", SEQ_ORDINAL), data_id = str_c(STATE, "-", SEQ_ORDINAL), fill = STATE), data = dataEvents) +

    ###
    geom_vline(colour="#aaaaaa", xintercept = arrangeRank$vlinePos1-0.5, size=0.3, alpha=0.9) +
    geom_vline(colour="#aaaaaa", xintercept = arrangeRank$vlinePos2-0.5, size=0.3, alpha=0.9) +
    geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    annotate("text", x=arrangeRank$vlinePos1-2, y=-1, label=arrangeRank$vlinePos1) +
    annotate("text", x=arrangeRank$vlinePos2-2, y=-1, label=arrangeRank$vlinePos2) +
    geom_hline(colour="#aaaaaa", yintercept = validFirstRank-1, size=0.3, alpha=0.9) +
    geom_hline(colour="#aaaaaa", yintercept = validLastRank, size=0.3, alpha=0.9) +
    annotate("text", y=middleValid, x=min(dataEvents$REL_STATE_START_DATE)-6, label=sprintf("N=%i \n  within \n timeframe" , nValid )) + # (%.2f%%) \n , percValid

    theme_classic() +
    scale_fill_manual(values=colorsDef) +
    labs(title = sprintf("%s N=%i (%.2f%% of starting with: %s.. N=%i (%.2f%% of the whole data set)).
                         \n Sorted by distance from: %s N=%i (%.2f %% of selected)
                         \n %i subjects have disatance less than %i from (%s) event %s to event %s
                         ",
                         selectedDrug,
                         length(unique(alignTime$SUBJECT_ID)),
                         length(unique(alignTime$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100,

                         pathClicked,
                         length(unique(selectedPaths$SUBJECT_ID)),
                         length(unique(selectedPaths$SUBJECT_ID))/length(unique(patPaths$SUBJECT_ID))*100,

                         arrBy,
                         length(unique(dataEvents$SUBJECT_ID)),
                         length(unique(dataEvents$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100,

                         nValid,
                         vline,
                         befAft,
                         selectedDrug,
                         arrBy
    )) +

      # geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    # geom_vline(colour="#EE5E00", xintercept = vline, size=0.3, alpha=0.9) +
    # annotate("text", x=vline-2, y = -0.5, label=vline) +
    # theme_classic() +
    # scale_fill_manual(values=colorsDef) +
    # labs(title = sprintf("Selected start: %s - ... N=%i (%.2f%% of the whole data set)
    #                      \n Aligned by: %s N=%i (%.2f%% of selected)
    #                      \n Sorted by distance from: %s N=%i (%.2f %% of selected)",
    #
    #                      pathClicked,
    #                      length(unique(selectedPaths$SUBJECT_ID)),
    #                      length(unique(selectedPaths$SUBJECT_ID))/length(unique(patPaths$SUBJECT_ID))*100,
    #
    #                      selectedDrug,
    #                      length(unique(alignTime$SUBJECT_ID)),
    #                      length(unique(alignTime$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100,
    #
    #                      arrBy,
    #                      length(unique(dataEvents$SUBJECT_ID)),
    #                      length(unique(dataEvents$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100)) +

    theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="black"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title=element_text(size=10)) +
    scale_x_continuous(name=sprintf("Time (days) before or after the STATE %s", selectedDrug)) +
    scale_y_continuous(name=" ")+
    guides(fill = guide_legend(
      title = "STATE",
      title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
      label.theme = element_text(face = "plain", size = 9, colour = "darkgrey")
    ))

  girafe(ggobj = p1,
         height_svg = 6,
         width_svg = 8,
         options = list(
           opts_tooltip(
             opacity = .8,
             css = "background-color:yellow; color:black; padding:2px; border-radius:5px;"
           ),
           opts_toolbar(position = "topright"),
           opts_selection(type = "single", only_shiny = FALSE, css = "stroke:#fdfd66;stroke-width:1.4"),
           opts_hover(css = "fill:#444444;stroke:#444444;cursor:pointer;")
         ))
}

#
# id_new = trajectoryDataPrep2(inputData)
# id2_new = trajectoryDataPrep2(inputData2)
#
# id_old = trajectoryDataPrep(inputData)
#
#
# patStateLevel = id_new[[3]]
# pathClicked = "PAP---Colposcopy"
# patPaths = id_new[[2]]
# colorsDef = id_new[[4]]
# selectedDrug = "HSIL-1"
# arrBy = "PAP-1"
# vline = 1000
#
# drugLevel02(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug, arrBy, vline)


#' Function plotting funnel of selections made - how many patients get filtered out from initial
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug

funnel <- function(patStateLevel, pathClicked, patPaths, selectedDrug, arrBy, vline, befAft){ #Alignment plot WITH out of cohorts

  selectedPaths <- patPaths %>%
    filter(str_detect(path1, str_c("^", pathClicked)))

  selectedLevels <- patStateLevel %>%
    filter(SUBJECT_ID %in% selectedPaths$SUBJECT_ID) %>%
    group_by(SUBJECT_ID) %>%
    mutate(REL_STATE_START_DATE = as.numeric(STATE_START_DATE - min(STATE_START_DATE))) %>%
    mutate(REL_STATE_END_DATE = as.numeric(STATE_END_DATE - min(STATE_START_DATE))) %>%
    ungroup()

  alignTime = selectedLevels %>%
    group_by(SUBJECT_ID) %>%
    filter(str_c(STATE, "-", SEQ_ORDINAL) == selectedDrug) %>%
    summarise(pathAlign = min(REL_STATE_START_DATE)) %>%
    ungroup()

  arrangeRank = selectedLevels %>%
    group_by(SUBJECT_ID) %>%
    filter(str_c(STATE, "-", SEQ_ORDINAL) == arrBy) %>%
    summarise(pathStart = min(REL_STATE_START_DATE)) %>%
    ungroup() %>%
    inner_join(alignTime, by = join_by(SUBJECT_ID)) %>%
    mutate(Rank1 = rank(pathAlign - pathStart, ties.method = "first"))

  arrangeRank <- arrangeRank %>%
    mutate(distance = pathStart-pathAlign)%>%
    filter(!if_any(everything(), is.na))


  if (befAft == "Before"){
    arrangeRank$vlinePos1 <- -vline
    arrangeRank$vlinePos2 <- 0
  } else if (befAft == "After"){
    arrangeRank$vlinePos1 <- 0
    arrangeRank$vlinePos2 <- vline
  } else {
    arrangeRank$vlinePos1 <- -vline
    arrangeRank$vlinePos2 <- vline
  }
  #arrangeRank$vlinePos <- vline

  validSubjects <- arrangeRank %>%
    filter(distance >= vlinePos1 & distance <= vlinePos2)
  nValid <- length(validSubjects$SUBJECT_ID)
  if (nValid == 0){
    validFirstRank <- 0
    validLastRank <- 0
  } else {
    validFirstRank <- min(validSubjects$Rank1)
    validLastRank <- max(validSubjects$Rank1)
  }

  dataEvents <- inner_join(arrangeRank, selectedLevels, by = join_by(SUBJECT_ID)) %>%
    mutate(REL_STATE_START_DATE = REL_STATE_START_DATE - pathAlign) %>%
    mutate(REL_STATE_END_DATE = REL_STATE_END_DATE - pathAlign)

  all_IDs <- length(unique(patPaths$SUBJECT_ID))
  selected_IDs <- length(unique(selectedPaths$SUBJECT_ID))
  align_IDs <- length(unique(alignTime$SUBJECT_ID))
  dist_IDs <- length(unique(dataEvents$SUBJECT_ID))

  fig <- plot_ly()
  fig <- fig %>%
    add_trace(type = "funnel",
              y = c(sprintf("1. Dataset size:"),
                    sprintf("2. Trajectory starts with %s:", pathClicked),
                    sprintf("3. In this group %s is present:", selectedDrug),
                    sprintf("4. In this group %s is present:", arrBy),
                    sprintf("5. Disatance less than %i (between %s and %s):", vline, selectedDrug, arrBy)),
              x = c(all_IDs, selected_IDs, align_IDs, dist_IDs, nValid),
              textposition = "inside",
              textinfo = "number of IDs in group",
              opacity = 0.65,
              marker = list(color = c("#0dc5c1", "#00a0c1", "#0080c1", "#0060c1", "#003088")),
              connector = list(line = list(color = "darkgrey", width = 3)))
  fig
}


#' Function plotting patient trajectories with OUT OF COHORT periods and aligning the trajectories by selected drug
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug
# alignArrangeAll <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug, arrBy, vline){ #Alignment plot WITH out of cohorts
#
#   ### added 5th April 2023
#   pathArrange1 <- patStateLevel %>%
#     select(SUBJECT_ID, TIME, DrugLEVEL) %>%
#     group_by(SUBJECT_ID) %>%
#     filter(DrugLEVEL == selectedDrug) %>%
#     filter(TIME == first(TIME)) %>%
#     mutate(pathAlign = TIME) %>%
#     select(SUBJECT_ID, pathAlign)
#
#   pathArrange2 <- patStateLevel %>%
#     select(SUBJECT_ID, TIME, DrugLEVEL) %>%
#     group_by(SUBJECT_ID) %>%
#     filter(DrugLEVEL == arrBy) %>%
#     filter(TIME == first(TIME)) %>%
#     mutate(pathStart = TIME) %>%
#     select(SUBJECT_ID, pathStart) %>%
#     left_join(pathArrange1)
#
#   pathArrange2 <- pathArrange2 %>%
#     mutate(distance = pathAlign - pathStart)
#   pathArrange2$distance[is.na(pathArrange2$distance)] <- min(pathArrange2$distance)-1 #0
#
#   pathArrange2 <- pathArrange2 %>%
#     select(SUBJECT_ID, distance) %>%
#     arrange(desc(distance))
#   pathArrange2$Rank1 <- 1:nrow(pathArrange2)
#   patStateLevel <- left_join(patStateLevel, pathArrange2)
#
#   selectedDrug <- selectedDrug
#   data <- patStateLevel
#   dataDrugs <- data %>%
#     group_by(SUBJECT_ID) %>%
#     mutate(MONTH1 = 1:n()) %>%
#     ungroup()
#   dataT1 <- dataDrugs %>%
#     filter(DrugLEVEL==selectedDrug) %>%
#     group_by(SUBJECT_ID) %>%
#     mutate(T1 = min(MONTH1)) %>%
#     ungroup() %>%
#     select(SUBJECT_ID, T1) %>%
#     unique()
#   dataDrugs <- left_join(dataDrugs, dataT1)
#   dataDrugs$TX <- dataDrugs$MONTH1 - dataDrugs$T1
#
#   # Set geom parameters
#   col_name <- dataDrugs$TX
#   row_name <- dataDrugs$Rank1 ## reordering change
#   colors <- colorsDef
#
#   p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
#     geom_tile_interactive() +
#     geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
#     geom_vline(colour="#EE5E00", xintercept = vline, size=0.3, alpha=0.9) +
#     annotate("text", x=vline-2, y=-1, label=vline) +
#     theme_classic() +
#     scale_fill_manual(values=colors) +
#     labs(title = sprintf("Align and sort with same conditions from the whole data set")) +
#     theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
#           axis.ticks.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.text.y = element_text(face="plain", size=8, colour="black"),
#           axis.ticks.y = element_blank(),
#           axis.line.y = element_blank(),
#           plot.title=element_text(size=10)) +
#     theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
#           axis.ticks.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.text.y = element_text(face="plain", size=8, colour="black"),
#           axis.ticks.y = element_blank(),
#           axis.line.y = element_blank()) +
#     scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
#     scale_y_continuous(name=" ")+
#     guides(fill = guide_legend(
#       title = "STATE",
#       title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
#       label.theme = element_text(face = "plain", size = 9, colour = "darkgrey")
#     ))
#
#   girafe(ggobj = p1,
#          height_svg = 6,
#          width_svg = 8,
#          options = list(
#            opts_tooltip(
#              opacity = .8,
#              css = "background-color:yellow; color:black; padding:2px; border-radius:5px;"
#            ),
#            opts_toolbar(position = "topright"),
#            opts_selection(type = "single", only_shiny = FALSE, css = "stroke:#fdfd66;stroke-width:1.4"),
#            opts_hover(css = "fill:#444444;stroke:#444444;cursor:pointer;")
#
#          ))
# }
#
