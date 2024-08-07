% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.R
\name{weight_rollmean_sub}
\alias{weight_rollmean_sub}
\title{Subfunction of weight_rollmean}
\usage{
weight_rollmean_sub(mat, k)
}
\arguments{
\item{mat}{matrix of normalized/adjusted counts}

\item{k}{size of window for weighting}
}
\value{
A matrix of smoothed counts
}
\description{
Take in a counts matrix and apply pyramidal weighting
with a window size k to create smoothed expression intensities
}
