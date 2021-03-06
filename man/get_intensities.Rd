% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load-and-filter.R
\name{get_intensities}
\alias{get_intensities}
\title{Get intensities}
\usage{
get_intensities(pg, sample_name = NULL)
}
\arguments{
\item{pg}{ProteinGroups data.frame.}
}
\description{
Extract the intensity values and put slice numbers as the column names.
}
\examples{
proteinGroups_path <-
system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
pg <- load_MQ(proteinGroups_path)
ints <- get_intensities(pg)
}
