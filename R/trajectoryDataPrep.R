#' Function preparing labels, colors and intermediate data frames "patStateLevel", "patPaths" and "freqPaths"
#'
#' @param inputData data read in when running function "trajectoryViz"
trajectoryDataPrep <- function(inputData){
  ### Find labels and create color palette  ----
  colnames(inputData)[which(names(inputData) == "STATE_LABEL")]<- "STATE"
  inputData$STATE <- gsub('-','',inputData$STATE)
  labels <- c(unique(inputData$STATE), "End") # "End" is required for RSunburst package input
  labels <- labels[!labels %in% c("START", "EXIT", "OUT OF COHORT", "End")]

  colors_all <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73",
                  "#0072B2", "#D55E00", "#CC79A7", "#004949",
                  "#b66dff", "#924900", "#ff6db6", "#490092", "#601A4A" )

  colors <- colors_all[1:length(labels)]
  colorsDef <- setNames(c(colors,"#cccccc"), c(labels, "OUT OF COHORT"))

  ## Data prep for "patStateLevel" as an input to plot patient-level pathways ----
  patStateLevel <- inputData %>% ## Label ALL LEVELS
    filter(STATE != "START" & STATE != "EXIT") %>%
    group_by(SUBJECT_ID) %>%
    arrange(SUBJECT_ID, STATE_START_DATE) %>%
    mutate(TIME = 1:n()) %>%
    ungroup() %>%
    group_by(SUBJECT_ID) %>%
    mutate(OtherLEVEL = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
    ungroup() %>%
    select(-STATE_END_DATE) %>%
    select(-STATE_START_DATE) %>%
    mutate(STATEOtherLEVEL = paste(STATE, OtherLEVEL, sep = "-"))
  data1 <- patStateLevel %>% ### Label DRUG LEVELS (without OUT OF COHORT)
    filter(STATE != "OUT OF COHORT") %>%
    group_by(SUBJECT_ID) %>%
    mutate(LEVEL = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
    ungroup()%>%
    mutate(STATELEVEL = paste(STATE, LEVEL, sep = "-"))
  patStateLevel <- left_join(patStateLevel, data1)
  patStateLevel1 <- patStateLevel %>%
    select(SUBJECT_ID, STATE, STATEOtherLEVEL) %>%
    unique() %>%
    arrange(STATE) %>%
    group_by(SUBJECT_ID, grp = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
    mutate(counts = 1:n()) %>%
    ungroup() %>%
    arrange(SUBJECT_ID) %>%
    select(-grp, -STATE)
  patStateLevel <- left_join(patStateLevel, patStateLevel1)
  patStateLevel <- patStateLevel %>%
    mutate(DrugLEVEL = paste(STATE, counts, sep = "-")) %>%
    select(SUBJECT_ID, TIME, STATE, OtherLEVEL, STATEOtherLEVEL, LEVEL, STATELEVEL, DrugLEVEL)

  #patStateLevel --> Table LEVELS in Shiny
  ## SUBJECT_ID + STATE +
  ## TIME(count from data set entrance, 1 piece of time period according to the data)
  ## OtherLEVEL(each consecutive block of a state ia a Level and Out Of cohort is considered as a level too)
  ## LEVEL (OUT of Cohort is removed from counting levels)

  ### Data prep for "patPaths" as an input to plot patient-level pathways   ----
  patStateLens<- patStateLevel %>% # Get Level Lengths to add to patPaths
    select(SUBJECT_ID, STATE, TIME, DrugLEVEL) %>%
    arrange(SUBJECT_ID) %>%
    filter(STATE != "OUT OF COHORT") %>%
    group_by(SUBJECT_ID, grp = with(rle(DrugLEVEL), rep(seq_along(lengths), lengths))) %>%
    mutate(COUNT_CONS = 1:n()) %>%
    ungroup() %>%
    select(-grp) %>%
    group_by(SUBJECT_ID, DrugLEVEL) %>%
    mutate(eraLen = max(COUNT_CONS)) %>%
    ungroup() %>%
    select(-TIME, -COUNT_CONS) %>%
    distinct() %>%
    group_by(SUBJECT_ID) %>%
    summarise(Len = paste0(eraLen, collapse = "-")) %>%
    ungroup() %>%
    separate(Len, c("Len1", "Len2", "Len3", "Len4", "Len5", "Len6", "Len7", "Len8", "Len9", "Len10", "Len11"), "-", convert = TRUE, fill="right") %>%
    replace_na(list(Len1=0, Len2=0, Len3=0, Len4=0, Len5=0, Len6=0, Len7=0, Len8=0, Len9=0, Len10=0, Len11=0))

  ### Data prep for "freq" - path FREQUENCES (sunburst input ) ----
  patPaths <- patStateLevel %>%
    filter(STATE != "OUT OF COHORT") %>%
    group_by(SUBJECT_ID, grp = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
    mutate(COUNT_CONS = 1:n()) %>%
    ungroup() %>%
    select(-grp) %>%
    filter(COUNT_CONS == 1) %>%
    select(-TIME) %>%
    select(-OtherLEVEL) %>%
    select(-COUNT_CONS) %>%
    group_by(SUBJECT_ID) %>%
    summarise(path = paste0(STATE, collapse = "-")) %>%
    ungroup() %>%
    group_by(SUBJECT_ID) %>%
    mutate(path = paste0(path, "-End", collapse = "-")) %>%
    ungroup() %>%
    mutate(path1 = gsub("-", " ", path)) %>%
    mutate(pathSep = path) %>%
    separate(pathSep, c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11"), "-", convert = TRUE, fill="right")
  patPaths <- left_join(patPaths, patStateLens)

  patPaths2 <- patStateLevel %>% # needed to add paths "OUT OF COHORT-End"
    filter((STATE != "OUT OF COHORT" | (STATE == "OUT OF COHORT" & OtherLEVEL == 1))) %>%
    group_by(SUBJECT_ID, grp = with(rle(STATE), rep(seq_along(lengths), lengths))) %>%
    mutate(COUNT_CONS = 1:n()) %>% ## find start of drug
    ungroup() %>%
    select(-grp) %>%
    filter(COUNT_CONS == 1) %>%
    select(-TIME) %>%
    select(-OtherLEVEL) %>%
    select(-COUNT_CONS) %>%
    group_by(SUBJECT_ID) %>%
    summarise(path = paste0(STATE, collapse = "-")) %>%
    ungroup() %>%
    group_by(SUBJECT_ID) %>%
    mutate(path = paste0(path, "-End", collapse = "-")) %>%
    ungroup() %>%
    mutate(path1 = gsub("-", " ", path)) %>%
    mutate(pathSep = path) %>%
    separate(pathSep, c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11"), "-", convert = TRUE, fill="right")
  patPaths2 <- left_join(patPaths2, patStateLens) %>%
    filter(path == "OUT OF COHORT-End")
  patPaths <- rbind(patPaths, patPaths2) %>%
    arrange(SUBJECT_ID)

  ## Count path frequences
  freqPaths <- patPaths %>%
    group_by(path) %>%
    summarise(freq = n()) %>%
    ungroup() %>%
    arrange(desc(freq))

  freqPaths2Percent <- freqPaths
  freqPaths2Percent <- freqPaths2Percent %>%
    mutate(percentage = round((freq/length(patPaths$SUBJECT_ID)*100),2))
  ### Return ----
  return(list(freqPaths, patPaths, patStateLevel, colorsDef, labels, colors, freqPaths2Percent))
}
