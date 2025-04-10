% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msm_frydman.R
\name{msm_frydman}
\alias{msm_frydman}
\title{Determine NPMLE for Multi State illness death Markov model using Frydman (1995)}
\usage{
msm_frydman(data, tol = 1e-08)
}
\arguments{
\item{data}{A \code{data.frame} containing the columns named:
\describe{
  \item{\code{delta}:}{Did a transition from 1 -> 2 occur? (binary: 0 = no, 1 = yes); 
  In the left-truncated case, delta = 2 indicates initially observed in state 2.}
  \item{\code{Delta}:}{Was the transition to state 3 observed? (binary: 0 = no, 1 = yes);}
  \item{\code{L}:}{Left timepoint of interval censored transition to state 2 (numeric);}
  \item{\code{R}:}{Right timepoint of interval censored transition to state 2 (numeric);}
  \item{\code{time}:}{Time of event (transition to 3) or right-censoring in state 2 (numeric);}
  \item{\code{trunc}:}{(optional) Left-truncation time (numeric); Only used for entries with delta = 2.}
}}

\item{tol}{Tolerance of the EM algorithm. Algorithm will stop when the absolute difference 
between current mass estimates and new estimates is smaller than the tolerance}
}
\value{
A list with the following entries:
\describe{
  \item{\code{data_idx}: }{A list containing the data used for the fit (\code{matdata}),
  the indices for which group a subject belongs to (\code{GroupX_idx}), some computational 
  parameters (see Frydman(1995)) and the unique failure times of the 2->3 and 1->3
  transitions respectively in \code{t_n_star} and \code{e_k_star};}
  \item{\code{supportMSM}: }{A list containing all transition intervals in \code{A} and
  the theoretical support intervals in \code{Q_mat};}
  \item{\code{z_lambda}: }{Computational quantities, see Frydman(1995);}
  \item{\code{cdf}: }{A list of functions that allow to recover the cdf for 
  the 1->3 (\code{F13}) and 1->2 (\code{F12}) transition and the cumulative hazard 
  for the 2->3 (\code{Lambda23}) transition.;}
}
}
\description{
Determine NPMLE for Multi State illness death Markov model using Frydman (1995)
}
\details{
For an illness death model (1 = healthy, 2 = ill, 3 = dead) estimate 
the NPMLE in the following form:
\describe{
  \item{\code{F12}:}{Cumulative distribution function of 1->2 transition;}
  \item{\code{F13}:}{Cumulative distribution function of 1->3 transition;}
  \item{\code{Lambda23}:}{Cumulative intensity of 2->3 transition;}
}
}
\examples{
data <- data.frame(delta = c(0, 0, 1, 1), Delta = c(0, 1, 0, 1),
                   L = c(NA, NA, 1, 1.5), R = c(NA, 3, 2, 3),
                   time = c(4, 5, 6, 7))

mod_frydman <- msm_frydman(data)
visualise_data(data, mod_frydman)


}
\references{
Frydman, H. (1995). Nonparametric Estimation of a Markov 
'Illness-Death' Process from Interval- Censored Observations, with 
Application to Diabetes Survival Data. Biometrika, 82(4), 773-789. 
\doi{10.2307/2337344}
}
