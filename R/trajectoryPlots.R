#' Function plotting patient trajectories with OUT OF COHORT periods
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
clustPlotsO <- function(patStateLevel, pathClicked, patPaths, colorsDef) {
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
  col_name <- data$TIME
  colors <- colorsDef

  plot <- ggplot(data, aes(x = col_name, y = Rank1, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    theme_classic() +
    scale_fill_manual(values = colors) +
    labs(title = sprintf("Selected: %s ... N = %i ", pathClicked, length(selected$SUBJECT_ID))) +
    theme(
      axis.text.x = element_text(face = "plain", size = 8, colour = "grey"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(face = "plain", size = 8, colour = "grey"),
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
         height_svg = 8,
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
      axis.text.x = element_text(face = "plain", size = 8, colour = "grey"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(face = "plain", size = 8, colour = "grey"),
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
         height_svg = 8,
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
drugLevel0 <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug){
  selectedDrug <- selectedDrug
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
  row_name <- dataDrugs$Rank1
  colors <- colorsDef

  p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    theme_classic() +
    scale_fill_manual(values=colors) +
    labs(title = sprintf("%s ... (%s: N= %i) ", pathClicked, selectedDrug, length(unique(dataT1$SUBJECT_ID)))) +
    theme(axis.text.x = element_text(face="plain", size=8, colour="grey"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="grey"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
    scale_y_continuous(name=" ")+
    guides(fill = guide_legend(
      title = "STATE",
      title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
      label.theme = element_text(face = "plain", size = 8, colour = "darkgrey")
    ))

  girafe(ggobj = p1,
         height_svg = 8,
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

#' Function plotting patient trajectories without OUT OF COHORT periods and aligning the trajectories by selected drug
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
#' @param selectedDrug reactive input from clustPlots plot for aligning the trajectories by selected drug
drugLevel <- function(patStateLevel, pathClicked, patPaths, colorsDef, selectedDrug){
  selectedDrug <- selectedDrug
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
  row_name <- dataDrugs$Rank1
  colors <- colorsDef

  p1 <- ggplot(dataDrugs, aes(x = col_name, y = row_name, tooltip = DrugLEVEL, data_id = DrugLEVEL, fill = STATE)) +
    geom_tile_interactive() +
    geom_vline(colour="#444444", xintercept = -0.5, size=0.6, alpha=0.9) +
    theme_classic() +
    scale_fill_manual(values=colors) +
    theme(axis.text.x = element_text(face="plain", size=8, colour="grey"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(face="plain", size=8, colour="grey"),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    labs(title = sprintf("%s ... (%s: N= %i) ", pathClicked, selectedDrug, length(unique(dataT1$SUBJECT_ID)))) +
    scale_x_continuous(name=sprintf("Times before or after the STATE %s", selectedDrug)) +
    scale_y_continuous(name="")+
    guides(fill = guide_legend(
      title = "STATE",
      title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
      label.theme = element_text(face = "plain", size = 8, colour = "darkgrey")
    ))

  girafe(ggobj = p1,
         height_svg = 8,
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

#' Function plotting patient trajectories with OUT OF COHORT periods arranged by length of the trajectory
#'
#' @param patStateLevel data frame from trajectoryDataPrep
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame from trajectoryDataPrep
#' @param colorsDef from trajectoryDataPrep
selectPatientsStartAlpha <- function(patStateLevel, pathClicked, patPaths, colorsDef) {
  data <- patStateLevel
  selected <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  data1 <- data %>%
    filter(SUBJECT_ID %in% selected$SUBJECT_ID)
  pat <- data1 %>%
    group_by(SUBJECT_ID) %>%
    summarise(sample_size = n()) %>%
    ungroup() %>%
    arrange(desc(sample_size)) %>%
    mutate(Rank = row_number()) %>%
    select(-sample_size)
  pat
  length(unique(pat$SUBJECT_ID))
  data <- left_join(data1, pat)
  # Set geom parameters
  col_name <- data$TIME
  row_name <- data$Rank
  colors <- colorsDef
  ggplot(data, aes(x = col_name, y = row_name, fill = STATE)) +
    geom_tile(aes(height = 1)) +
    theme_classic() +
    scale_fill_manual(values = colors) +
    theme(
      axis.text.x = element_text(face = "plain", size = 8, colour = "grey"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(face = "plain", size = 8, colour = "grey"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    scale_x_continuous("Time in the data set") +
    scale_y_continuous("Patients") +
    guides(fill = guide_legend(
      title = "STATE",
      title.theme = element_text(size = 10, face = "plain", colour = "darkgrey", angle = 0),
      label.theme = element_text(face = "plain", size = 8, colour = "darkgrey")
    ))
}

