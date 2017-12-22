% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetworkToolbox--master.R
\name{convertConnBrainMat}
\alias{convertConnBrainMat}
\title{Import CONN Toolbox Brain Matrices to R format and convert to correlations}
\usage{
convertConnBrainMat(MatlabData = file.choose(), progBar = TRUE)
}
\arguments{
\item{MatlabData}{Input for Matlab data file. Defaults to interactive file choice}

\item{progBar}{Should progress bar be displayed? Defaults to TRUE. Set FALSE for no progress bar}
}
\value{
Returns an array of correlation connectivity (n x n x m)
}
\description{
Converts a Matlab brain z-score connectivity file (n x n x m) where \strong{n} is the n x n connectivity matrices and \strong{m} is the participant
}
\examples{
\dontrun{braindata<-convertConnBrainMat()}
}
\author{
Alexander Christensen <alexpaulchristensen@gmail.com>
}