% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectoryPlots.R
\name{drugLevel0}
\alias{drugLevel0}
\title{Function plotting patient trajectories with OUT OF COHORT periods and aligning the trajectories by selected drug}
\usage{
drugLevel0(
  patStateLevel,
  pathClicked,
  patPaths,
  colorsDef,
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

\item{colorsDef}{from trajectoryDataPrep}

\item{selectedDrug}{reactive input from clustPlots plot for aligning the trajectories by selected drug}

\item{arrBy}{reactive user input}

\item{vline}{reactive user input}
}
\description{
Function plotting patient trajectories with OUT OF COHORT periods and aligning the trajectories by selected drug
}
