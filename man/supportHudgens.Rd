% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supportHudgens.R
\name{supportHudgens}
\alias{supportHudgens}
\title{Determine the support of the NPMLE for interval censored data.}
\usage{
supportHudgens(intervals, reduction = TRUE, existence = FALSE)
}
\arguments{
\item{intervals}{A data.frame with 3 columns containing half-open intervals (left open, right closed)
and an indicator whether the interval results from a censored transition
 or truncation:
\describe{
  \item{\code{L}:}{Left side of interval;}
  \item{\code{R}:}{Right side of interval;}
  \item{\code{cens}:}{Indicator whether interval resulted from interval censoring or left truncation
  (1 = censoring, 0 = truncation);}
  \item{\code{id}:}{(optional) Identifier for the observation this interval belongs to (numeric/integer).
   Only required if existence = TRUE;}
} Note that the truncation intervals need to be in the form (N, Inf] with N a numeric value.}

\item{reduction}{Should the support be reduced using Lemma 3 from Hudgens (2005)? This 
requires checking an extra condition. Default is TRUE.}

\item{existence}{Should the existence of the NPMLE be checked using Theorem 1/Lemma 4 from 
Hudgens (2005)? Requires \code{id} to be present in \code{intervals}. Default is FALSE.}
}
\value{
\itemize{
\item \code{graph}: An \code{igraph} object representing the censoring/truncation intervals
\item \code{support}: Support estimated from the censoring intervals
\item \code{dir_graph}: A directed \code{igraph} object used to determine whether the NPMLE
exists in the presence of left-truncation.
\item \code{exist_mle}: Logical output indicating whether the NPMLE exists.
}
}
\description{
Given censoring/truncation intervals, find the maxcliques and determine
 the support of the interval censored problem.
}
\references{
Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
 \doi{10.1111/j.1467-9868.2005.00516.x}
}
