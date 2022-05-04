#' Average length or times of 1st STATE level
#'
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame with all patients and their paths
stats1 <- function(pathClicked, patPaths) {
  selected <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  sprintf("%s: %.2f (%.0f) ", selected$L1[1], mean(as.integer(selected$Len1)), median(as.integer(selected$Len1)))
}
#' Average length or times of 2nd STATE level
#'
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame with all patients and their paths
stats2 <- function(pathClicked, patPaths) {
  selected <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  sprintf("%s: %.2f (%.0f) ", selected$L2[1], mean(as.integer(selected$Len2)), median(as.integer(selected$Len2)))
}
#' Average length or times of 3rd STATE level
#'
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame with all patients and their paths
stats3 <- function(pathClicked, patPaths) {
  selected <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  sprintf("%s: %.2f (%.0f) ", selected$L3[1], mean(as.integer(selected$Len3)), median(as.integer(selected$Len3)))
}
#' Average length or times of 4th STATE level
#'
#' @param pathClicked reactive input from clicking on the sunburst
#' @param patPaths data frame with all patients and their paths
stats4 <- function(pathClicked, patPaths) {
  selected <- patPaths %>%
    filter(str_sub(path1, 1, nchar(pathClicked)) == pathClicked)
  sprintf("%s: %.2f (%.0f) ", selected$L4[1], mean(as.integer(selected$Len4)), median(as.integer(selected$Len4)))
}

