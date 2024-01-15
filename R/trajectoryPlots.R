#' Function cutting data from selected start
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param startState reactive input from user input
cutStart <- function(patStateLevel, startState){
  selectedStart <- patStateLevel %>%
    select(SUBJECT_ID, STATE, TIME, DrugLEVEL)
  selectedStart1 <- selectedStart %>%
    filter(DrugLEVEL == startState) %>% # selectedDrug + -1
    mutate(selectedStartTIME = selectedStart$TIME[selectedStart$DrugLEVEL == startState]) %>%
    select(SUBJECT_ID, selectedStartTIME) %>%
    unique() %>%
    group_by(SUBJECT_ID) %>%
    filter(selectedStartTIME == min(selectedStartTIME))
  #n_distinct(selectedStart1$SUBJECT_ID)

  selectedStart <- left_join(selectedStart, selectedStart1)
  #n_distinct(selectedStart$SUBJECT_ID)

  selectedStartData <- selectedStart %>%
    filter(selectedStart$TIME >= selectedStart$selectedStartTIME) %>%
    group_by(SUBJECT_ID) %>%
    arrange(SUBJECT_ID, TIME) %>%
    mutate(STATE_START_DATE = 1:n()) %>%
    mutate(STATE_END_DATE = 1:n()) %>%
    ungroup()#%>%
    #select(-TIME, -DrugLEVEL, -selectedStartTIME)

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
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)

  mat_data <- selected %>%
    select(SUBJECT_ID, Len1, Len2, Len3)
  mat <- mat_data %>%
    as.data.frame() %>%
    column_to_rownames("SUBJECT_ID") %>%
    as.matrix()

  hc <- hclust(dist(mat, method = "euclidean"), method = "complete")
  ord1 <- tibble(
    Rank1 = 1:nrow(mat),
    SUBJECT_ID = as.double(rownames(mat)[hc$order])
  )

  data1 <- data %>%
    filter(SUBJECT_ID %in% selected$SUBJECT_ID)
  data <- left_join(data1, ord1)

  dataDrugs <- data %>%
    filter(STATE != "OUT OF COHORT") %>%
    group_by(SUBJECT_ID) %>%
    mutate(MONTH1 = 1:n()) %>%
    ungroup()
  dataDrugs <- left_join(dataDrugs, patPaths)

  # Set geom parameters
  col_name1 <- dataDrugs$MONTH1
  colors <- colorsDef
  plot <- ggplot(dataDrugs, aes(x = col_name1, y = Rank1, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    theme_classic() +
    scale_fill_manual(values = colors) +
    labs(title = sprintf("Selected: %s ... N = %i ", pathClicked, length(selected$SUBJECT_ID))) +
    theme(
      axis.text.x = element_text(face = "plain", size = 8, colour = "black"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(face = "plain", size = 8, colour = "black"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous("Time (unit period depends on input data)") +
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

#' Function plotting patient trajectories with OUT OF COHORT periods and aligning the trajectories by selected drug
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug
#' @param arrBy reactive user input
#' @param vline reactive user input

drugLevel0 <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug, arrBy, vline){ #Alignment plot WITH out of cohorts

  selectedPaths <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  selectedLevels <- patStateLevel %>%
    filter(SUBJECT_ID %in% selectedPaths$SUBJECT_ID)

  alignTime <- selectedLevels %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == selectedDrug) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathAlign = TIME) %>%
    select(SUBJECT_ID, pathAlign)

  arrangeRank <- selectedLevels %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == arrBy) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathStart = TIME) %>%
    select(SUBJECT_ID, pathStart) %>%
    left_join(alignTime)

  arrangeRank <- arrangeRank %>%
    mutate(distance = pathAlign - pathStart)%>%
    filter(!if_any(everything(), is.na)) %>%
    select(SUBJECT_ID, distance, pathAlign) %>%
    arrange(desc(distance))
  arrangeRank$Rank1 <- 1:nrow(arrangeRank)
  dataDrugs <- left_join(arrangeRank, selectedLevels)
  dataDrugs$TX <- dataDrugs$TIME - dataDrugs$pathAlign

  # Set geom parameters
  col_name <- dataDrugs$TX
  row_name <- dataDrugs$Rank1
  colors <- colorsDef

  p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    geom_vline(colour="#EE5E00", xintercept = vline, size=0.3, alpha=0.9) +
    annotate("text", x=vline-2, y=-1, label=vline) +
    theme_classic() +
    scale_fill_manual(values=colors) +
    labs(title = sprintf("Selected start: %s - ... N=%i (%.2f%% of the whole data set)
                         \n Aligned by: %s N=%i (%.2f%% of selected)
                         \n Sorted by distance from: %s N=%i (%.2f %% of selected)",

                         pathClicked,
                         length(unique(selectedPaths$SUBJECT_ID)),
                         length(unique(selectedPaths$SUBJECT_ID))/length(unique(patPaths$SUBJECT_ID))*100,

                         selectedDrug,
                         length(unique(alignTime$SUBJECT_ID)),
                         length(unique(alignTime$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100,

                         arrBy,
                         length(unique(dataDrugs$SUBJECT_ID)),
                         length(unique(dataDrugs$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100)) +

    theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="black"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title=element_text(size=10)) +
    scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
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


#' Function plotting funnel of selections made - how many patients get filtered out from initial
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug

funnel <- function(patStateLevel, pathClicked, patPaths, selectedDrug, arrBy){ #Alignment plot WITH out of cohorts

  selectedPaths <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  selectedLevels <- patStateLevel %>%
    filter(SUBJECT_ID %in% selectedPaths$SUBJECT_ID)

  alignTime <- selectedLevels %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == selectedDrug) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathAlign = TIME) %>%
    select(SUBJECT_ID, pathAlign)

  arrangeRank <- selectedLevels %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == arrBy) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathStart = TIME) %>%
    select(SUBJECT_ID, pathStart) %>%
    left_join(alignTime)

  arrangeRank <- arrangeRank %>%
    mutate(distance = pathAlign - pathStart)%>%
    filter(!if_any(everything(), is.na)) %>%
    select(SUBJECT_ID, distance, pathAlign) %>%
    arrange(desc(distance))
  arrangeRank$Rank1 <- 1:nrow(arrangeRank)
  dataDrugs <- left_join(arrangeRank, selectedLevels)
  #ataDrugs$TX <- dataDrugs$TIME - dataDrugs$pathAlign

  # Set geom parameters
  #col_name <- dataDrugs$TX
  #row_name <- dataDrugs$Rank1
  #colors <- colorsDef

  # p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
  #   geom_tile_interactive() +
  #   geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
  #   geom_vline(colour="#EE5E00", xintercept = vline, size=0.3, alpha=0.9) +
  #   annotate("text", x=vline-2, y=-1, label=vline) +
  #   theme_classic() +
  #   scale_fill_manual(values=colors) +
  #   labs(title = sprintf("Selected start: %s - ... N=%i (%.2f%% of the whole data set)
  #                        \n Aligned by: %s N=%i (%.2f%% of selected)
  #                        \n Sorted by distance from: %s N=%i (%.2f %% of selected)",
  #
  #                        pathClicked,
  #                        length(unique(selectedPaths$SUBJECT_ID)),
  #                        length(unique(selectedPaths$SUBJECT_ID))/length(unique(patPaths$SUBJECT_ID))*100,
  #
  #                        selectedDrug,
  #                        length(unique(alignTime$SUBJECT_ID)),
  #                        length(unique(alignTime$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100,
  #
  #                        arrBy,
  #                        length(unique(dataDrugs$SUBJECT_ID)),
  #                        length(unique(dataDrugs$SUBJECT_ID))/length(unique(selectedPaths$SUBJECT_ID))*100)) +
  #
  #   theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
  #         axis.ticks.x = element_blank(),
  #         axis.line.x = element_blank(),
  #         axis.text.y = element_text(face="plain", size=8, colour="black"),
  #         axis.ticks.y = element_blank(),
  #         axis.line.y = element_blank(),
  #         plot.title=element_text(size=10)) +
  #   scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
  #   scale_y_continuous(name=" ")+
  #   guides(fill = guide_legend(
  #     title = "STATE",
  #     title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
  #     label.theme = element_text(face = "plain", size = 9, colour = "darkgrey")
  #   ))
  #
  # girafe(ggobj = p1,
  #        height_svg = 6,
  #        width_svg = 8,
  #        options = list(
  #          opts_tooltip(
  #            opacity = .8,
  #            css = "background-color:yellow; color:black; padding:2px; border-radius:5px;"
  #          ),
  #          opts_toolbar(position = "topright"),
  #          opts_selection(type = "single", only_shiny = FALSE, css = "stroke:#fdfd66;stroke-width:1.4"),
  #          opts_hover(css = "fill:#444444;stroke:#444444;cursor:pointer;")
  #        ))

  all_IDs <- length(unique(patPaths$SUBJECT_ID))
  selected_IDs <- length(unique(selectedPaths$SUBJECT_ID))
  align_IDs <- length(unique(alignTime$SUBJECT_ID))
  dist_IDs <- length(unique(dataDrugs$SUBJECT_ID))

  fig <- plot_ly()
  fig <- fig %>%
    add_trace(type = "funnel",
              y = c("1. All Subjects", "2. Selected from Sunburst", "3. Selected by Aligning", "4. Have the State of Distance"),
              x = c(all_IDs, selected_IDs, align_IDs, dist_IDs),
              textposition = "inside",
              textinfo = "number of IDs in group",
              opacity = 0.65,
              marker = list(color = c("#0dc5c1", "#00a0c1", "#0080c1", "#0060c1")),
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
alignArrangeAll <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug, arrBy, vline){ #Alignment plot WITH out of cohorts

  ### added 5th April 2023
  pathArrange1 <- patStateLevel %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == selectedDrug) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathAlign = TIME) %>%
    select(SUBJECT_ID, pathAlign)

  pathArrange2 <- patStateLevel %>%
    select(SUBJECT_ID, TIME, DrugLEVEL) %>%
    group_by(SUBJECT_ID) %>%
    filter(DrugLEVEL == arrBy) %>%
    filter(TIME == first(TIME)) %>%
    mutate(pathStart = TIME) %>%
    select(SUBJECT_ID, pathStart) %>%
    left_join(pathArrange1)

  pathArrange2 <- pathArrange2 %>%
    mutate(distance = pathAlign - pathStart)
  pathArrange2$distance[is.na(pathArrange2$distance)] <- min(pathArrange2$distance)-1 #0

  pathArrange2 <- pathArrange2 %>%
    select(SUBJECT_ID, distance) %>%
    arrange(desc(distance))
  pathArrange2$Rank1 <- 1:nrow(pathArrange2)
  patStateLevel <- left_join(patStateLevel, pathArrange2)

  selectedDrug <- selectedDrug
  data <- patStateLevel
  dataDrugs <- data %>%
    group_by(SUBJECT_ID) %>%
    mutate(MONTH1 = 1:n()) %>%
    ungroup()
  dataT1 <- dataDrugs %>%
    filter(DrugLEVEL==selectedDrug) %>%
    group_by(SUBJECT_ID) %>%
    mutate(T1 = min(MONTH1)) %>%
    ungroup() %>%
    select(SUBJECT_ID, T1) %>%
    unique()
  dataDrugs <- left_join(dataDrugs, dataT1)
  dataDrugs$TX <- dataDrugs$MONTH1 - dataDrugs$T1

  # Set geom parameters
  col_name <- dataDrugs$TX
  row_name <- dataDrugs$Rank1 ## reordering change
  colors <- colorsDef

  p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    geom_vline(colour="#EE5E00", xintercept = vline, size=0.3, alpha=0.9) +
    annotate("text", x=vline-2, y=-1, label=vline) +
    theme_classic() +
    scale_fill_manual(values=colors) +
    labs(title = sprintf("Align and sort with same conditions from the whole data set")) +
    theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="black"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title=element_text(size=10)) +
    theme(axis.text.x = element_text(face="plain", size=8, colour="black"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="black"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
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

