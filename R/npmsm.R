#' NPMLE for general multi-state model with interval censored transitions 
#' 
#' @description For a general Markov chain multi-state model with interval censored 
#' transitions calculate the NPMLE of the transition intensities. The estimates 
#' are returned as an \code{\link[mstate:msfit]{msfit}} object.
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
#' \code{c("multinomial", "poisson")}, with multinomial the default. Multinomial 
#' will use the EM algorithm described in Gomon and Putter (2024) and Poisson
#' will use the EM algorithm described in Gu et al. (2023).
#' @param inits Which distribution should be used to generate the initial estimates 
#' of the intensities in the EM algorithm. One of c("equalprob", "unif", "beta"), 
#' with "equalprob" assigning 1/K to each intensity, with K the number of distinct 
#' observation times (\code{length(unique(gd[, "time"]))}). For "unif", each 
#' intensity is sampled from the Unif[0,1] 
#' distribution and for "beta" each intensity is sampled from the Beta(a, b) distribution. 
#' If "beta" is chosen, the argument \code{beta_params} must be specified as a 
#' vector of length 2 containing the parameters of the beta distribution.
#' Default = "equalprob".
#' @param beta_params A vector of length 2 specifying the beta distribution parameters 
#' for initial distribution generation. First entry will be used as \code{shape1}
#' and second entry as \code{shape2}. See \code{help(rbeta)}. Only used if \code{inits = "beta"}. 
#' @param support_manual Used for specifying a manual support region for the transitions.
#' A list of length the number of transitions in \code{tmat}, 
#' each list element containing a data frame with 2 named columns L and R indicating the 
#' left and right values of the support intervals. When specified, all intensities 
#' outside of these intervals will be set to zero for the corresponding transitions.
#' Intensities set to zero cannot be changed by the EM algorithm. Will use inits = "equalprob".
#' @param exact Numeric vector indicating to which states transitions are observed at exact times.
#' Must coincide with the column number in \code{tmat}.
#' @param maxit Maximum number of iterations. Default = 100.
#' @param tol Tolerance of the convergence procedure. A change in the value of 
#' \code{conv_crit} in an iteration of less than \code{tol} will make the procedure stop.
#' @param conv_crit Convergence criterion. Stops procedure when the difference 
#' in the chosen quantity between two consecutive iterations is smaller 
#' than the tolerance level \code{tol}. One of the following:
#' \describe{
#' \item{"haz"}{Stop when change in maximum estimated intensities (hazards) \code{ < tol}.}
#' \item{"prob"}{Stop when change in estimated probabilities \code{ < tol}.}
#' \item{"lik"}{Stop when change in observed-data likelihood \code{ < tol}.}
#' } Default is "haz". The options "haz" and "lik" can be compared across different
#' \code{method}s, but "prob" is dependent on the chosen \code{method}. Most 
#' conservative (requiring most iterations) is "prob", followed by "haz" and finally "lik".
#' @param verbose Should iteration messages be printed? Default is FALSE
#' @param manual Manually specify starting transition intensities? If \code{TRUE},
#' the transition intensity for each bin for each transition must be entered manually.
#' DO NOT USE for large data sets, and in general it is not adviced to use this.
#' @param newmet Should contributions after last observation time also be used 
#' in the likelihood? Default is FALSE.
#' @param include_inf Should an additional bin from the largest observed time to 
#' infinity be included in the algorithm? Default is FALSE.
#' @param checkMLE Should a check be performed whether the estimate has converged 
#' towards a true Maximum Likelihood Estimate? This is done by comparing 
#' the reduced gradient to the value of \code{checkMLE_tol}. Default is TRUE.
#' @param checkMLE_tol Tolerance for checking whether the estimate has converged to MLE.
#' Whenever an estimated transition intensity is smaller than \code{prob_tol}, it is assumed 
#' to be zero and its reduced gradient is not considered for determining whether 
#' the NPMLE has been reached. Default = \code{1e-4}.
#' @param prob_tol If an estimated probability is smaller than \code{prob_tol}, 
#' it will be set to zero during estimation. Default value is \code{tol/10}.
#' @param remove_redundant Should redundant observations be removed before running 
#' the algorithm? An observation is redundant when the same state has been observed 
#' more than 3 times consecutively, or if it is a repeat observation of an 
#' absorbing state. Default is TRUE.
#' @param remove_bins Should a bin be removed during the algorithm if all
#' estimated intensities are zero for a single bin? Can improve 
#' computation speed for large data sets. Note that zero means the estimated intensities 
#' are smaller than \code{prob_tol}. Default is FALSE.
#' @param estimateSupport Should the support of the transitions be estimated using 
#' the result of Hudgens (2005)? Currently produces incorrect support sets - 
#' DO NOT USE. Default = \code{FALSE}
#' @param init_int A vector of length 2, with the first entry indicating what 
#' percentage of mass should be distributed over (second entry) what percentage 
#' of all first bins. Default is c(0, 0), in which case the argument is ignored.
#' This argument has no practical uses and only exists for demonstration purposes 
#' in the related article.
#' @param ... Further arguments to \code{\link{estimate_support_msm}}
#' 
#' 
#' @return A list with the following entries:
#' \describe{
#'   \item{\code{A}: }{A list of class \code{\link[mstate:msfit]{msfit}} containing 
#'   the cumulative intensities for each transition and the transition matrix used;}
#'   \item{\code{Ainit}: }{Initial intensities, in an object of class \code{\link[mstate:msfit]{msfit}};}
#'   \item{\code{gd}: }{Data used for the estimation procedure;}
#'   \item{\code{ll}: }{Log-likelihood value of the procedure at the last iteration;}
#'   \item{\code{delta}: }{Change in log-likelihood value at the last iteration;}
#'   \item{\code{it}: }{Number of iterations of the procedure;}
#'   \item{\code{taus}: }{Unique time points of the data, the cumulative intensity
#'   only makes jumps at these time points.;}
#'   \item{\code{tmat}: }{The transition matrix used, see \code{\link[mstate:transMat]{transMat}};}
#'   \item{\code{tmat2}: }{A summary of the transitions in the model, see \code{\link[mstate:to.trans2]{to.trans2}};}
#'   \item{\code{ll_history}: }{The log-likelihood value at every iteration of the procedure;}
#'   \item{\code{KKT_violated}: }{How often were KKT conditions violated during
#'   maximisation of the likelihood? In other words, how often did we hit the optimization 
#'   boundary during the procedure?;}
#'   \item{\code{min_time}: }{The smallest time of an observation in the used data. 
#'   Note that the smallest time in the data is used as zero reference;}
#'   \item{\code{reduced_gradient}: }{The reduced gradient at the last iteration.
#'   Rows indicate the transitions and columns the unique observation times;}
#'   \item{\code{isMLE}: }{Has the procedure converged to the NPMLE? Checked 
#'   using \code{checkMLE_tol};}
#'   \item{\code{langrangemultiplier}: }{The lagrange multipliers at the last iteration;}
#'   \item{\code{aghmat}: }{A matrix representation of the transition intensities in \code{A}.
#'   Rows represent transitions and columns unique observation times;}
#'   \item{\code{Ygk}: }{The summed at-risk indicator for all subjects in the data at the last iteration.
#'   Rows represent transitions and columns unique observation times;}
#'   \item{\code{Dmk}: }{The summed probability of making a transition for all subjects at the last iteration.
#'   Rows represent transitions and columns unique observation times;}
#'   \item{\code{method}: }{Method used for the optimization procedure;}
#'   \item{\code{maxit}: }{Maximum number of allowed iterations;}
#'   \item{\code{tol}: }{Tolerance of the convergence procedure;}
#'   \item{\code{conv_crit}: }{Convergence criterion of the procedure;}
#'   \item{\code{checkMLE}: }{Was the reduced gradient checked at the last iteration to determine convergence?;}
#'   \item{\code{checkMLE_tol}: }{The tolerance of the checkMLE procedure;}
#'   \item{\code{prob_tol}: }{Tolerance for probabilities to be set to zero;}
#'   \item{\code{remove_redundant}: }{Were redundant observations removed before performing the procedure?;}
#' }
#' 
#' 
#' 
#' @details
#' Denote the unique observation times in the data as \eqn{0 = \tau_0, \tau_1, \ldots, \tau_K}{0 = tau_0, tau_1, ..., tau_K}
#' Let \eqn{g, h \in H}{g, h in H} denote the possible states in the model and \eqn{X(t)}{X(t)} the state of the process at time t.
#' 
#' Then this function can be used to estimate the transition intensities 
#' \eqn{\alpha_{gh}^k = \alpha_{gh}(\tau_k)}{alpha_(gh)^k = alpha_(gh)(tau_k)}. 
#' 
#' Having obtained these estimated, it is possible to recover the transition probabilities 
#' \eqn{\mathbf{P}(X(t) = h | X(s) = g)}{P(X(t) = h | X(s) = g)} for \eqn{t > s}{t > s} using 
#' the \code{\link{transprob}} functions.
#' 
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @importFrom igraph is_dag graph_from_adjacency_matrix
#' @import checkmate
#' @export
#' 
#' @seealso \code{\link{transprob}} for calculating transition probabilities,
#'  \code{\link{plot.npmsm}} for plotting the cumulative intensities, 
#'  \code{\link{print.npmsm}} for printing some output summaries,
#'  \code{\link{visualise_msm}} for visualising data,
#'  \code{\link[mstate:msfit]{msfit}} for details on the output object.
#' 
#' 
#' @references   
#' D. Gomon and H. Putter, Non-parametric estimation of transition intensities 
#' in interval censored Markov multi-state models without loops,
#' arXiv, 2024, \doi{10.48550/arXiv.2409.07176}
#' 
#' Y. Gu, D. Zeng, G. Heiss, and D. Y. Lin, 
#' Maximum likelihood estimation for semiparametric regression models with 
#' interval-censored multistate data, Biometrika, Nov. 2023, \doi{10.1093/biomet/asad073}
#' 
#' Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#' 
#' 
#' @examples 
#' #Create transition matrix using mstate functionality: illness-death
#' if(require(mstate)){
#'   tmat <- mstate::trans.illdeath()
#' }
#' 
#' #Write a function for evaluation times: observe at 0 and uniform inter-observation times.
#' eval_times <- function(n_obs, stop_time){
#'   cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
#' }
#' 
#' #Use built_in function to simulate illness-death data
#' #from Weibull distributions for each transition
#' sim_dat <- sim_id_weib(n = 50, n_obs = 6, stop_time = 15, eval_times = eval_times,
#' start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
#' 
#' tmat <- mstate::trans.illdeath()
#' 
#' #Fit the model using method = "multinomial"
#' mod_fit <- npmsm(gd = sim_dat, tmat = tmat, tol = 1e-2)
#' 
#' #Plot the cumulative intensities for each transition
#' plot(mod_fit)
#' 



npmsm <- function(gd, tmat, method = c("multinomial", "poisson"),
                  inits = c("equalprob", "homogeneous", "unif", "beta"), beta_params, support_manual, 
                  exact, maxit = 100, tol = 1e-4, conv_crit = c("haz", "prob", "lik"),
                  verbose = FALSE, manual = FALSE, 
                  newmet = FALSE, include_inf = FALSE, checkMLE = TRUE,
                  checkMLE_tol = 1e-4, prob_tol = tol/10,
                  remove_redundant = TRUE, remove_bins = FALSE,
                  estimateSupport = FALSE, init_int = c(0, 0), ...){
  

# Argument Checks ---------------------------------------------------------

  
  method <- match.arg(method, choices = c("multinomial", "poisson"))
  conv_crit <- match.arg(conv_crit, choices = c("haz", "prob", "lik"))
  inits <- match.arg(inits, choices = c("equalprob", "homogeneous", "unif", "beta"))
  
  
  arg_checks <- makeAssertCollection()
  assertMatrix(tmat, add = arg_checks)
  tmat2 <- mstate::to.trans2(tmat)
  assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  assertNames(names(gd), must.include = c("id", "state", "time"), 
              add = arg_checks)
  assertIntegerish(gd[, "id"], add = arg_checks)
  assertIntegerish(gd[, "state"], add = arg_checks)
  assertNumeric(gd[, "time"], lower = 0, add = arg_checks)
  if(!missing(exact)){
    assertIntegerish(exact, lower = 1, upper = nrow(tmat), add = arg_checks)
  }
  if(!missing(support_manual)){
    assertList(support_manual, types = c("matrix", "list"), len = nrow(tmat2), add = arg_checks)
  }
  if(inits == "beta"){
    assertNumeric(beta_params, lower = 0, upper = Inf, any.missing = FALSE, len = 2, add = arg_checks)
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
    out <- EM_multinomial( gd = gd, tmat = tmat, tmat2 = tmat2, inits = inits, beta_params = beta_params, support_manual = support_manual,
                        exact = exact, maxit = maxit, tol = tol, conv_crit = conv_crit, manual = manual, 
                        verbose = verbose, newmet = newmet, include_inf = include_inf, checkMLE = checkMLE,
                        checkMLE_tol = checkMLE_tol, prob_tol = prob_tol, remove_bins = remove_bins, init_int = init_int,
                        ... )  
  } else if( method == "poisson" ){
    out <- EM_poisson( gd = gd, tmat = tmat, tmat2 = tmat2, inits = inits, beta_params = beta_params, support_manual = support_manual,
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