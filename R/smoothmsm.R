#' Smooth hazard estimation for general multi-state model with interval censored
#' transitions 
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
#' @param inits Which distribution should be used to generate the initial estimates 
#' of the intensities in the EM algorithm. One of c("equalprob", "unif", "beta"), 
#' with "equalprob" assigning 1/K to each intensity, with K the number of distinct 
#' observation times (\code{length(unique(gd[, "time"]))}). For "unif", each 
#' intensity is sampled from the Unif[0,1] 
#' distribution and for "beta" each intensity is sampled from the Beta(a, b) distribution. 
#' If "beta" is chosen, the argument \code{beta_params} must be specified as a 
#' vector of length 2 containing the parameters of the beta distribution.
#' Default = "equalprob".
#' @param exact Numeric vector indicating to which states transitions are observed at exact times.
#' Must coincide with the column number in \code{tmat}.
#' @param formula Formula to interpret in data for covariates.
#' @param data A \code{data.frame} containing a column called \code{'id'} 
#' (identifying the subjects in \code{gd}) and 
#' variables which to interpret in \code{formula}.
#' @param deg_splines Degree to use for the B-spline basis functions. Defaults 
#' to 3 (cubic B-splines).
#' @param n_segments Number of segments to use for the P-splines. The
#' segments will space the domain evenly. According to Eilers \& Marx (2021), this 
#' number cannot be chosen too large. Default = 20. 
#' @param ord_penalty Order of the P-spline penalty (penalty on the difference 
#' between d-order differences of spline coefficients). See Eilers \& Marx 
#' Section 2.3. Defaults to 2.
#' @param n_bins How many evenly spaced bins should the domain be split into 
#' for the estimation problem (unrelated to the splines). The evenly spaced bins 
#' must not contain any two consecutive transition intervals of a single subject.
#' By default, the number of bins is chosen equal to the length of the 
#' smallest observed transition interval.
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
#' @param prob_tol If an estimated probability is smaller than \code{prob_tol}, 
#' it will be set to zero during estimation. Default value is \code{tol/10}.
#' 
#' 
#' 
#' @references 
#' Eilers, P.H.C. and Marx, B.D., Practical Smoothing: The Joys of P-splines, 
#' Cambridge University Press (2021)
#' 
#' 
#' 
#' @importFrom mstate to.trans2 msfit probtrans
#' @importFrom igraph is_dag graph_from_adjacency_matrix
#' @import checkmate
#' @export
#' 
#' 
#' 
#' 
#' 



smoothmsm <- function(gd, tmat, exact, formula, data,
                      deg_splines = 3, n_segments = 20, ord_penalty = 2, n_bins,
                      maxit = 100, tol = 1e-4, conv_crit = c("haz", "prob", "lik"),
                      verbose = FALSE, prob_tol = tol/10){
  
  

# Pre-processing ---------------------------------------


  ## Argument checks ---------------------------------------------------------
  call <- match.call()
  conv_crit <- match.arg(conv_crit, choices = c("haz", "prob", "lik"))
  
  arg_checks <- makeAssertCollection()
  assertMatrix(tmat, add = arg_checks)
  assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  assertNames(names(gd), must.include = c("id", "state", "time"), 
              add = arg_checks)
  assert_class(gd[, "id"], classes = c("numeric", "character"), add = arg_checks)
  assertIntegerish(gd[, "state"], add = arg_checks)
  assertNumeric(gd[, "time"], lower = 0, add = arg_checks)
  if(!missing(exact)){
    assertIntegerish(exact, lower = 1, upper = nrow(tmat), add = arg_checks)
  }
  assertIntegerish(maxit, lower = 0, add = arg_checks)
  assertNumeric(tol, add = arg_checks)
  assertLogical(verbose, add = arg_checks)
  
  #Check if all states in gd are also in tmat
  assertTRUE(all(unique(gd[, "state"]) %in% 1:nrow(tmat2)), add = arg_checks)
  
  # Check Spline arguments
  
  #B-spline basis must be at least degree 0
  assertIntegerish(deg_splines, lower = 0, upper = Inf, max.len = 1, min.len = 1,
                   add = arg_checks)
  #P-spline penalty order must be at least 1
  assertIntegerish(ord_penalty, lower = 1, upper = Inf, min.len = 1, max.len = 1,
                   add = arg_checks)
  
  assertFormula(formula, null.ok = TRUE, add = arg_checks)
  

  ## Data processing ---------------------------------------------------------

  #Data transformations for probtrans compatibility
  #We set the smallest time to 0, because probtrans always calculates from 0 onward.
  min_time <- min(gd[, "time"])
  gd[which(gd$time == min_time), "time"] <- 0
  
  max_time <- max(gd[, "time"])
  
  tmat2 <- mstate::to.trans2(tmat)
  
  #TODO: NEED TO CHECK WHETHER the given n_intervals can be used
  #Must be larger than the number if we were to use interval size equal to smallest transition interval
  #Note that domain_length = max_time, as we set smallest time to 0
  
  ## Derive parameters from data ---------------------------------------------
  
  subject_names <- unique(gd[, "id"])
  n_subjects <- length(subject_names)

  n_states <- nrow(tmat) # no of states
  n_transitions <- nrow(tmat2) # no of transitions

  # Spline quantities
  n_splines <- deg_splines + n_segments
  
  

  ## Derived Quantities Check ------------------------------------------------

  #If covariates specified.
  if(!is.null(data)){
    subject_names_data <- unique(data[, "id"])
    
    #Check whether all observed subjects also have related covariates
    assertTRUE(all(subject_names %in% subject_names_data), add = arg_checks)
  }
  
  #See Eilers/Marx (2021) Section 2.3
  assertTRUE(ord_penalty < n_splines, add = arg_checks)
  
  #Need to perform checks on spline derivatives
  #Order must be smaller than number of B-spline basis functions
  #Number of B-spline basis functions is equal to degree of B-splines + number of segments
  #Seeing as segments are evenly spaced
  
  
  #Report assertions
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  
  
  ## Additional parameter determination ------------------------------------
  
  penalization_matrix <- diff(diag(n_splines), diff=ord_penalty)

# EM Algorithm ------------------------------------------------------------

  
  #Pass by reference using environments:
  #https://bookdown.org/content/d1e53ac9-28ce-472f-bc2c-f499f18264a3/reference.html
  
  #Initiate initial estimates
  
  #Initiate 


  ## E step ------------------------------------------------------------------



  ## M Step ------------------------------------------------------------------

  
  
    
  
  
  
  
}