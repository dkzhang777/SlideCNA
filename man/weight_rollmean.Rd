% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{weight_rollmean}
\alias{weight_rollmean}
\title{Expressional smoothing along a chromosome using a weighted pyramidal moving average}
\usage{
weight_rollmean(dat, k = 101)
}
\arguments{
\item{dat}{data.table of normalized/adjusted counts}

\item{k}{size of window for weighting}
}
\value{
A data.table of expression intensities
}
\description{
Take in a data.table of genomic positions and bead normalized/modified counts and apply pyramidal weighting
with a window size k to create smoothed expression intensities
}
