% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support_from_direct_intervals.R
\name{support_from_direct_intervals}
\alias{support_from_direct_intervals}
\title{Estimate support of multiple transitions given direct transition intervals}
\usage{
support_from_direct_intervals(direct_intervals, tmat)
}
\arguments{
\item{direct_intervals}{Output from \code{\link{direct_from_observed_intervals}}.}

\item{tmat}{A transition matrix as created by \code{transMat}}
}
\value{
A list containing the estimated support sets for each possible transition
in \code{tmat}.
}
\description{
Given only direct transition intervals, determine the support
for each transition separately using Hudgens(2001) result. Each state is
considered from a competing-risks viewpoint. Hudgens(2005) result is applied 
to see if the NPMLE for any of the transitions does not exist.
}
