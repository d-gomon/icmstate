% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probtrans_coxph.R
\name{expand_covariates_long_data}
\alias{expand_covariates_long_data}
\title{Expand covariates for a data frame so that covariates can be transition 
specific.}
\usage{
expand_covariates_long_data(newdata)
}
\arguments{
\item{newdata}{A \code{data.frame} containing a single row for each transition per
subject in the data. For a model with m possible transitions, and n subjects 
\code{newdata} must have m*n rows. It must contain the following named columns:
\describe{
  \item{\code{id}:}{Unique identifier of the subject, can be numeric or character;}
  \item{\code{from}:}{State from which the transition takes place;}
  \item{\code{to}:}{State to which the transition takes place;}
  \item{\code{trans}:}{Transition number in the \code{'transMat'} \code{trans} this transition relates to;}
  \item{\code{"variables"}:}{The variables and their values for the subject
  identified by "id" for the transition this entry relates to. Names must 
  match the names of the variables in coxph \code{object}.}
} Note that newdata must contain a column containing the variable which was 
used to determine the stratum of a transition in \code{object}. 
Usually the stratum is determined from one of the required columns. The 
"variables" columns can usually be obtained using \code{\link[mstate:expand.covs]{expand.covs}}.}
}
\description{
Expand covariates for a data frame so that covariates can be transition 
specific.
}
\keyword{internal}
