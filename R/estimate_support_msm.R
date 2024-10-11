#' Estimate the support of a general Markov interval-censored Multi-state model without 
#' loops.
#' 
#' @description Given a realisation of a multi-state model, estimate the support 
#' of the different transitions possible in that MSM. The estimation is performed 
#' by viewing each possible state in a competing risks setting and applying the 
#' result of Hudgens (2001) to determine the support and left-truncation 
#' intervals and Hudgens (2005) to check whether a solution is possible.
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param tmat A transition matrix as created by \code{transMat}
#' 
#' @return TODO
#' 
#' 
#' @references Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#'
#' @references M. G. Hudgens, G. A. Satten, and I. M. Longini, 
#' Nonparametric Maximum Likelihood Estimation for Competing Risks Survival Data
#'  Subject to Interval Censoring and Truncation, Biometrics, vol. 57, no. 1, 
#'  Pages 74-80, March 2001, \doi{10.1111/j.0006-341x.2001.00074.x}
#' 
#' @import igraph
#' @importFrom mstate to.trans2
#' @export
#' 
#' 




estimate_support_msm <- function(gd, tmat){
  
  #This function operates using a few steps:
  #Step 1: Get all observed transitions from the data
  #Note that the observed transitions may contain indirect transitions, 
  #censored observations in a transient state and transitions which could 
  #have occurred in multiple ways. We therefore perform:
  #Step 2: For each observed transitions, determine all "worst-case" direct transitions
  #By doing so, we have assumed the least possible information from each
  #observed transition. This allows us to perform
  #Step 3: Apply Hudgens(2001) result to intervals from Step 2 to obtain an
  #estimate of the support for each possible direct transition.
  #Note that Hudgens(2005) can then be used to determine whether the NPMLE for
  #that transition exists or not. It is however unclear whether this result 
  #generalises to multi-state models, therefore should be used with caution.
  

# Get transition intervals from full data ---------------------------------
  #We obtain all transitions observed in the data in a data.frame with columns
  #entry_time   time_from   time_to   from   to   id
  observed_intervals <- get_trans_intervals(gd = gd, tmat = tmat)
  #Note that censored observations have 
  #entry_time: time first seen in censored state
  #time_from: last time seen in censored state
  #time_to: Inf
  

# Obtain full intervals into direct transition intervals --------
  
  #We want to obtain direct transition intervals from the full intervals
  #This means that any transition that either cannot have happened directly
  #or a transition that could have happened in different ways 
  #(e.g. the 1->3 transition in an illness-death model)
  #will be split up into multiple components
  
  direct_intervals <- direct_from_observed_intervals(observed_intervals = observed_intervals,
                                                                tmat = tmat, gd = gd)
    
  

# Determine support sets using Hudgens (2001) -----------------------------

  #From the direct transitions, we now want to determine subsets corresponding 
  #to each transition and feed each of these subsets with the correct
  #censoring and truncation intervals to supportHudgens()
  
  support_sets <- support_from_direct_intervals(direct_intervals = direct_intervals,
                                                tmat = tmat)
  
  
  return(support_sets)
}





