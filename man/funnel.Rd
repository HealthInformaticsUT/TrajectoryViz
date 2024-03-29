% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectoryPlots.R
\name{funnel}
\alias{funnel}
\title{Function plotting funnel of selections made - how many patients get filtered out from initial}
\usage{
funnel(
  patStateLevel,
  pathClicked,
  patPaths,
  selectedDrug,
  arrBy,
  vline,
  befAft
)
}
\arguments{
\item{patStateLevel}{data frame from trajectoryDataPrep}

\item{pathClicked}{reactive input from clicking on the sunburst}

\item{patPaths}{data frame from trajectoryDataPrep}

\item{selectedDrug}{reactive input from clustPlots plot for aligning the trajectories by selected drug}
}
\description{
Function plotting funnel of selections made - how many patients get filtered out from initial
}
