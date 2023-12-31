% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictions.R
\name{compute_shapley}
\alias{compute_shapley}
\title{Compute shapley values}
\usage{
compute_shapley(explainer, contribution_mat)
}
\arguments{
\item{explainer}{An \code{explain} object.}

\item{contribution_mat}{The contribution matrix.}
}
\value{
A \code{data.table} with shapley values for each test observation.
}
\description{
Compute shapley values
}
\keyword{internal}
