#' NPMLE for general multi-state model with interval censored transitions 
#' 
#' @description For a general Markov chain multi-state model with interval censored 
#' transitions calculate the NPMLE. 
#' 
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param tmat A transition matrix as created by \code{\link[mstate:transMat]{transMat}}
#' @param method Which method should be used for the EM algorithm. Choices are 
#' \code{c("multinomial", "poisson")}, with multinomial the default.
#' @param support_manual Used for specifying a manual support region for the transitions.
#' A list of length the number of transitions in \code{tmat}, 
#' each list element containing a data frame with 2 names columns L and R indicating the 
#' left and right values of the support intervals
#' @param exact Numeric vector indicating to which states transitions are observed at exact times.
#' Must coincide with the column number in \code{tmat}.
#' @param maxit Maximum number of iterations.
#' @param tol Tolerance of the procedure. A change in the value of 
#' \code{conv_crit} in an iteration of less than \code{tol} will make the procedure stop.
#' @param conv_crit Convergence criterion. Stops procedure when the difference 
#' in the chosen quantity between two consecutive iterations is smaller 
#' than the tolerance level \code{tol}. One of the following:
#' \describe{
#' \item{"haz"}{Stop when change in estimated intensities (hazards) \code{< tol}.}
#' \item{"prob"}{Stop when change in estimated probabilities \code{< tol}.}
#' \item{"lik"}{Stop when change in observed-data likelihood \code{< tol}.}
#' } Default is "haz". The options "haz" and "lik" can be compared across different
#' \code{method}s, but "prob" is dependent on the chosen \code{method}. Most 
#' conservative (safest) is "prob", followed by "haz" and finally "lik".
#' @param verbose Should iteration messages be printed? Default is FALSE
#' @param manual Manually specify starting transition intensities?
#' @param newmet Should contributions after last observation time also be used 
#' in the likelihood? Default is FALSE.
#' @param include_inf Should an additional bin from the largest observed time to 
#' infinity be included in the algorithm? Default is FALSE.
#' @param checkMLE Should a check be performed whether the estimate has converged 
#' towards a true Maximum Likelihood Estimate? Default is TRUE.
#' @param checkMLE_tol Tolerance for checking whether the estimate has converged to MLE.
#' Whenever an estimated transition intensity is smaller than the tolerance, it is assumed 
#' to be zero.
#' @param prob_tol If an estimated probability is smaller than \code{prob_tol}, 
#' it will be set to zero during estimation. Default value is \code{tol}.
#' @param remove_redundant Should redundant observations be removed before running 
#' the algorithm? Default is TRUE.
#' @param remove_bins Should bins be removed during the algorithm if there is 
#' 0 estimated intensity in all of the transitions? Significantly improves 
#' computation speed for large data sets. Note that 0 means the estimated intensities 
#' are smaller than \code{prob_tol}. Default is FALSE.
#' @param estimateSupport Should the support of the transitions be estimated using 
#' the result of Hudgens (2005)? Currently produces incorrect support sets - 
#' DO NOT USE.
#' @param init_int A vector of length 2, with the first entry indicating what 
#' percentage of mass should be distributed over (second entry) what percentage 
#' of all first bins. Default is c(0, 0), in which case the argument is ignored.
#' This argument has no practical uses and only exists for demonstration purposes.
#' @param ... Further arguments to \code{\link{estimate_support_msm}}
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @importFrom igraph is_dag
#' @import checkmate
#' @export
#' 
#' 
#' @references Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573–587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#' 
#' 
#' @examples 
#' #Create transition matrix using mstate functionality
#' if(require(mstate)){
#'   tmat <- mstate::trans.illdeath()
#' }
#' 
#' #Write a function for evaluation times
#' eval_times <- function(n_obs, stop_time){
#'   cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
#' }
#' 
#' #Use built_in function to simulate from Weibull distributions for each transition
#' sim_dat <- sim_id_weib(n = 50, n_obs = 6, stop_time = 15, eval_times = eval_times,
#' start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
#' 
#' tmat <- mstate::trans.illdeath()
#' mod_fit <- npmsm(gd = sim_dat, tmat = tmat, tol = 1e-2)
#' plot(mod_fit$A)
#' 



npmsm <- function(gd, tmat, method = c("multinomial", "poisson"), support_manual, 
                  exact, maxit = 100, tol = 1e-4, conv_crit = c("haz", "prob", "lik"),
                  verbose = FALSE, manual = FALSE, 
                  newmet = FALSE, include_inf = FALSE, checkMLE = TRUE,
                  checkMLE_tol = 1e-4, prob_tol = tol,
                  remove_redundant = TRUE, remove_bins = FALSE,
                  estimateSupport = FALSE, init_int = c(0, 0), ...){
  

# Argument Checks ---------------------------------------------------------

  
  method <- match.arg(method, choices = c("multinomial", "poisson"))
  conv_crit <- match.arg(conv_crit, choices = c("haz", "prob", "lik"))
  
  
  arg_checks <- makeAssertCollection()
  assertMatrix(tmat, add = arg_checks)
  tmat2 <- mstate::to.trans2(tmat)
  assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  assertNames(names(gd), must.include = c("id", "state", "time"), 
              add = arg_checks)
  assertIntegerish(gd[["id"]], add = arg_checks)
  assertIntegerish(gd[["state"]], add = arg_checks)
  assertNumeric(gd[["time"]], lower = 0, add = arg_checks)
  if(!missing(exact)){
    assertIntegerish(exact, lower = 1, upper = nrow(tmat), add = arg_checks)
  }
  if(!missing(support_manual)){
    assertList(support_manual, types = c("matrix", "list"), len = nrow(tmat2), add = arg_checks)
  }
  assertIntegerish(maxit, lower = 0, add = arg_checks)
  assertNumeric(tol, add = arg_checks)
  assertLogical(manual, add = arg_checks)
  assertLogical(verbose, add = arg_checks)
  
  
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  
  #Transitions in transmat (not their names) are always named by the corresponding column.
  #We make use of this a lot!
  

  #--------------- Check whether MSM contains loops ------------------#
  
  #We can check whether the MSM contains loops by checking whether the
  #graph corresponding to the adjacency matrix is acyclic (is a Directed Acyclic Graph (DAG)).
  contains_loops <- FALSE
  adjacency_matrix <- tmat
  adjacency_matrix[!is.na(tmat)] <- 1
  adjacency_graph <- graph_from_adjacency_matrix(adjacency_matrix)
  if(!is_dag(adjacency_graph)){
    contains_loops <- TRUE
    warning("Multi-state model contains cycles, NPMLE cannot be determined. 
            Algorithm will run, but results cannot be interpreted!")
  }
  
  
  

# Data Transformations ----------------------------------------------------
  
  #TO-DO: We need to check for loops in the state space, as when we have loops
  #the NPMLE might not be identifiable (it is unclear whether observing 
  #a transient state twice implies that we stayed in the state or that we 
  #transitioned back and forth)
  #To implement this: check for loops in remove_redundant_observations and 
  #throw warning if there are loops. If there are loops, do not remove any non-
  #absorbing observations.
  
  #Remove redundant observations to lighten the computational load
  if(isTRUE(remove_redundant)){
    gd <- remove_redundant_observations(gd = gd, tmat = tmat)  
  }
  
  #Make prob_tol at least 1e-4, otherwise some probabilities get set to zero too fast.
  prob_tol <- min(prob_tol, 1e-4)
  
  
  

# Determination of Support ------------------------------------------------
  #browser()
   if(estimateSupport){
     estimated_support <- estimate_support_msm(gd, tmat)
     support_manual <- lapply(estimated_support, '[[', 2)
     exist_mle <- all(sapply(estimated_support, '[[', 4))
   }

# Call to EM algorithm ----------------------------------------------------

  if( method == "multinomial" ){
    out <- EM_multinomial( gd = gd, tmat = tmat, tmat2 = tmat2, support_manual = support_manual,
                        exact = exact, maxit = maxit, tol = tol, conv_crit = conv_crit, manual = manual, 
                        verbose = verbose, newmet = newmet, include_inf = include_inf, checkMLE = checkMLE,
                        checkMLE_tol = checkMLE_tol, prob_tol = prob_tol, remove_bins = remove_bins, init_int = init_int,
                        ... )  
  } else if( method == "poisson" ){
    out <- EM_poisson( gd = gd, tmat = tmat, tmat2 = tmat2, support_manual = support_manual,
                       exact = exact, maxit = maxit, tol = tol, conv_crit = conv_crit, manual = manual, 
                       verbose = verbose, newmet = newmet, include_inf = include_inf, checkMLE = checkMLE,
                       checkMLE_tol = checkMLE_tol, prob_tol = prob_tol, remove_bins = remove_bins, init_int = init_int,
                       ... )
  }
  
  if(estimateSupport){
    out$estimated_support <- estimated_support
    out$exist_mle <- exist_mle
  }
  
  # Return additional run information
  out$method <- method
  if(!missing(exact)){
    out$exact <- exact  
  }
  out$maxit <- maxit
  out$tol <- tol
  out$conv_crit <- conv_crit
  out$checkMLE <- checkMLE
  out$checkMLE_tol <- checkMLE_tol
  out$prob_tol <- prob_tol
  out$remove_redundant <- remove_redundant
  
  
  class(out) <- "npmsm"
  return(out)
}