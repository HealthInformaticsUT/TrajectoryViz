harmonize_data = function(inputData){
  if("STATE_LABEL" %in% colnames(inputData)){
    inputData = inputData %>%
      rename("STATE" = "STATE_LABEL")
  }
  if("OUT OF COHORT" %in% unique(inputData$STATE)){
    inputData = inputData %>%
      group_by(SUBJECT_ID) %>%
      arrange(SUBJECT_ID, STATE_START_DATE) %>%
      group_by(SUBJECT_ID) %>%
      mutate(GROUP = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
      group_by(SUBJECT_ID, STATE, GROUP) %>%
      summarize(STATE_START_DATE = min(STATE_START_DATE), STATE_END_DATE = max(STATE_END_DATE)) %>%
      ungroup()
  }

  inputData = inputData %>%
    arrange(SUBJECT_ID, STATE_START_DATE) %>%
    group_by(SUBJECT_ID, STATE) %>%
    mutate(SEQ_ORDINAL = rank(STATE_START_DATE)) %>%
    ungroup() %>%
    filter(STATE != "START" & STATE != "EXIT" & STATE != "OUT OF COHORT") %>%
    select(SUBJECT_ID, STATE, STATE_START_DATE, STATE_END_DATE, SEQ_ORDINAL)

  inputData = inputData %>%
    mutate(STATE = str_replace_all(STATE, "-", " "))
  # Handle '+' cases in STATE_LABEL
  inputData = inputData %>%
    mutate(STATE = str_replace_all(STATE, "\\+", ""))

  return(inputData)
}

#' Function preparing labels, colors and intermediate data frames "patStateLevel", "patPaths" and "freqPaths"
#'
#' @param inputData data read in when running function "trajectoryViz"
trajectoryDataPrep <- function(inputData){
  patStateLevel = harmonize_data(inputData)

  # Patient stay lengths table
  patStateLens = patStateLevel %>%
    group_by(SUBJECT_ID) %>%
    mutate(eraLen = as.numeric(STATE_END_DATE - STATE_START_DATE)) %>%
    mutate(N = 1:n()) %>%
    ungroup() %>%
    pivot_wider(id_cols = SUBJECT_ID, names_from = N, names_prefix = "Len", values_from = eraLen, values_fill = 0)

  patPaths = patStateLevel %>%
    group_by(SUBJECT_ID) %>%
    mutate(N = 1:n()) %>%
    mutate(
      path = str_c(c(STATE, "End"), collapse = "-"),
      path1 = str_c(c(STATE, "End"), collapse = "---")
    ) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(SUBJECT_ID, path, path1), names_from = N, names_prefix = "L", values_from = STATE)

  patPaths = left_join(patPaths, patStateLens, by = "SUBJECT_ID")

  # Frequent paths
  freqPaths <- patPaths %>%
    group_by(path) %>%
    summarise(freq = n()) %>%
    ungroup() %>%
    arrange(desc(freq))

  freqPaths2Percent <- freqPaths
  freqPaths2Percent <- freqPaths2Percent %>%
    mutate(percentage = round((freq/length(patPaths$SUBJECT_ID)*100),2))

  # Colors & stuff
  labels <- c(unique(patStateLevel$STATE), "End") # "End" is required for RSunburst package input
  labels <- labels[!labels %in% c("START", "EXIT", "OUT OF COHORT", "End")]

  colors_all <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73",
                  "#0072B2", "#D55E00", "#CC79A7", "#004949",
                  "#b66dff", "#924900", "#ff6db6", "#490092", "#601A4A" )

  colors <- colors_all[1:length(labels)]
  colorsDef <- setNames(c(colors,"#cccccc"), c(labels, "OUT OF COHORT"))

  return(list(freqPaths, patPaths, patStateLevel, colorsDef, labels, colors, freqPaths2Percent))
}


#' Input data validation, sets columns to upper case and checks for mandatory columns
#'
#' @param inputData input to trajectoryViz without previous preprocessing
#'
#' @return list with feedback and input data with upper case column names
#' @internal
inputDataValidation <- function(inputData){
  feedback <- NULL
  feedback1 <- NULL
  feedback2 <- NULL
  data <- inputData
  
  # set colnames to upper
  colnames(data) <- toupper(colnames(inputData))
  
  # check for colnames 
  mandatoryCols <- c("SUBJECT_ID", "STATE_START_DATE", "STATE_END_DATE", "SEQ_ORDINAL") # "state_label"
  missingCols <- ""
  for (name in mandatoryCols){
    if (!(name %in% colnames(data))){
      missingCols <- paste(missingCols,name)
    } 
  }
  
  # check if state_label or state present 
  if((!("STATE" %in% colnames(data)) & !("STATE_LABEL" %in% colnames(data)))){
    missingCols <- paste(missingCols,"STATE_LABEL")
  }
  
  if (nchar(missingCols) > 0){
    feedback1 <- paste0("Missing mandatory columns:", missingCols)
  }
   
  # TODO add new line
  #subject_id
  if (typeof(data$SUBJECT_ID) != "integer" | typeof(data$SUBJECT_ID) != "double"){
    feedback2 <- paste0("Subject id column type is ", typeof(data$SUBJECT_ID), ", but must be integer or double. ")
  }
  #state_label
  if (typeof(data$STATE_LABEL) != "character"){
    feedback2 <- paste0(feedback2,"State label column type is ", typeof(data$STATE_LABEL), ", but must be character. ")
  }
  #state_start & state_end
  if (typeof(data$STATE_START_DATE) != "character"){
    feedback2 <- paste0(feedback2,"State start date column type is ", typeof(data$STATE_START_DATE), ", but must be character. ")
  }
  
  if (typeof(data$STATE_END_DATE) != "character"){
    feedback2 <- paste0(feedback2,"State end date column type is ", typeof(data$STATE_END_DATE), ", but must be character. ")
  }
  #seq_ordinal 
  if (typeof(data$STATE_START_DATE) != "integer"){
    feedback2 <- paste0(feedback2,"Seq ordinal date column type is ", typeof(data$SEQ_ORDINAL), ", but must be integer. ")
  }
  
  # feedback 
  feedback <- c(feedback1, feedback2)
  
  return(list(feedback,data))
  
}
