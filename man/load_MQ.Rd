% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load-and-filter.R
\name{load_MQ}
\alias{load_MQ}
\title{Load and filter MaxQuant proteinGroups.txt}
\usage{
load_MQ(proteinGroups_path, ignore_slices = NULL, sample_name = NULL)
}
\arguments{
\item{proteinGroups_path}{Path to proteinGroups.txt}
}
\description{
Load the results and remove contaminants.
}
\examples{
proteinGroups_path <-
system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
pg <- load_MQ(proteinGroups_path)
}
