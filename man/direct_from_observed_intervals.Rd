% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/direct_from_observed_intervals.R
\name{direct_from_observed_intervals}
\alias{direct_from_observed_intervals}
\title{Translate observed transition intervals into direct transition intervals}
\usage{
direct_from_observed_intervals(observed_intervals, tmat, gd)
}
\arguments{
\item{observed_intervals}{Output from \code{\link{get_trans_intervals}}.}

\item{tmat}{A transition matrix as created by \code{transMat}}

\item{gd}{A \code{data.frame} with the following named columns
\describe{
  \item{\code{id}:}{Subject idenitifier;}
  \item{\code{state}:}{State at which the subject is observed at \code{time};}
  \item{\code{time}:}{Time at which the subject is observed;}
} The true transition time between states is then interval censored between the times.}
}
\value{
A \code{data.frame} with the following named columns
\describe{
  \item{\code{entry_time}:}{Time of entry into "from" state;}
  \item{\code{time_from}:}{Last time subject(id) was seen in state "from";}
  \item{\code{time_to}:}{First time subject(id) was seen in state "to";}
  \item{\code{from}:}{State from which a transition was observed;}
  \item{\code{to}:}{State to which the transition was observed;}
  \item{\code{id}:}{Subject identifier;}
} For right-censored observations, entry_time denotes the first time
seen in the censored state, time_from the last time seen in the censored state,
time_to is \code{Inf}, from the censored state and to is \code{NA}.
}
\description{
Given observed transition intervals, determine the "worst" (least informative) possible 
direct transition intervals that could have occurred to form this sample.
}
