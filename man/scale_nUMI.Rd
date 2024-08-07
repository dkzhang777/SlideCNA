% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin.R
\name{scale_nUMI}
\alias{scale_nUMI}
\title{Scale for nUMI (UMI Count) to generate CNV scores}
\usage{
scale_nUMI(dat_bin, thresh_hard = FALSE)
}
\arguments{
\item{dat_bin}{data.table of relative expression intensities per bin}

\item{thresh_hard}{TRUE if using strict thresholds for expression thresholds and FALSE if adjusting
thresholds based on 1 + or - the mean of absolute min and max values}
}
\value{
data.table of CNV scores per bin
}
\description{
This function re-scales expression intensities to be in a smaller range, normalizes for nUMI per bin,
and subtracts reference bead signal
}
