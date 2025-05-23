% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_probtrans.R
\name{plot_probtrans}
\alias{plot_probtrans}
\title{Plot the transition probabilities for a fitted \code{npmsm} model}
\usage{
plot_probtrans(
  npmsmlist,
  from = NULL,
  to = NULL,
  transitions = NULL,
  landmark,
  interpolate = TRUE,
  facet = TRUE,
  times_interpol = NULL,
  c.legend = TRUE,
  c.names = NULL
)
}
\arguments{
\item{npmsmlist}{An "npmsm" object or a list containing multiple "npmsm" objects}

\item{from}{A numeric value indicating the state from which we consider the 
transition probabilities. Default is NULL, meaning we consider transition probabilities 
from all states from which a direct transition is possible.}

\item{to}{A numeric vector indicating to which states we consider the transition 
probabilities. Only states that can be reached from the \code{"from"} state are considered.}

\item{transitions}{A numeric vector indicating which transitions to consider (plot). 
Can only be used if \code{"from"} is not specified, as it only works for direct transitions.}

\item{landmark}{A landmark time indicating from which time on survival should be determined.
If missing, the smallest between the time in the first "npmsm" object or 0 will be used.}

\item{interpolate}{Should the cumulative hazard be linearly interpolated before 
determining transition probabilities? Default is TRUE.}

\item{facet}{Should the resulting plot be faceted (one panel per transition)? 
Default is TRUE.}

\item{times_interpol}{At which times should the cumulative hazard be interpolated?
Only necessary to specify if interpolate = TRUE.}

\item{c.legend}{Should legend be displayed for colour (different entries in 
\code{npmsmlist})? Default is TRUE.}

\item{c.names}{A character vector indicating the names to display in the legend.
These names should represent the entries in \code{npmsmlist.}. Default = \code{NULL}.}
}
\value{
A plot will be produced in the plotting window. When assigning  
the output to an object, the underlying data frame used for plotting 
and a \code{'ggplot'} object will be returned in a list.
}
\description{
For a fitted \code{npmsm} model plot the transition probabilities
from a certain state for all possible (direct and indirect) transitions.
}
\examples{
require(mstate)
require(ggplot2)
#Generate from an illness-death model with exponential transitions with 
#rates 1/2, 1/10 and 1 for 10 subjects over a time grid.
gd <- sim_weibmsm(tmat = trans.illdeath(), shape = c(1,1,1),
                  scale = c(2, 10, 1), n_subj = 10, obs_pars = c(2, 0.5, 20), 
                  startprobs = c(0.9, 0.1, 0))
#Fit 2 models: 1 with at most 4 iterations and 1 with at most 20
mod1 <- npmsm(gd, trans.illdeath(), maxit = 4)
mod2 <- npmsm(gd, trans.illdeath(), maxit = 20)

#Plot the transition probabilities from state 1, without interpolating 
#the cumulative hazard for the npmsm runs with max 4 and 20 iterations.
plot_probtrans(list(mod1, mod2), from = 1, interpolate = FALSE,
               c.names = c("4 iterations", "20 iterations"))


}
