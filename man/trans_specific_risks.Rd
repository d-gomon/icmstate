% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probtrans_coxph.R
\name{trans_specific_risks}
\alias{trans_specific_risks}
\title{Calculate subject specific risks for subjects in newdata}
\usage{
trans_specific_risks(object, newdata, trans)
}
\arguments{
\item{object}{A \code{\link[survival:coxph]{coxph}} object fit on multi-state
data. Must contain a \code{strata(X)} term. Data used for the coxph() fit preferably prepared 
using \code{\link[mstate:msprep]{msprep}}.}

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

\item{trans}{A transition matrix as created by \code{\link[mstate:transMat]{transMat}}.}
}
\description{
Return a 2 dimensional array, with subjects in the first dimension and
transition (numbers) in the second. Can expand later to include time-dependent
covariates by introducing extra dimension. Entries of the matrix are exp(lp)_i^m,
with i denoting the subject and m the transition.
}
\keyword{internal}
