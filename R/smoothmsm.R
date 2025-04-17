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
#' segments will space the domain evenly. According to Eilers \& Marx (2021), it 
#' it OK to choose this number very large. Default = 20. 
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

  if(missing(formula)){
    formula <- NULL
  }

  ## Argument checks ---------------------------------------------------------
  call <- match.call()
  conv_crit <- match.arg(conv_crit, choices = c("haz", "prob", "lik"))
  
  arg_checks <- makeAssertCollection()
  assertMatrix(tmat, add = arg_checks)
  assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  assertNames(names(gd), must.include = c("id", "state", "time"), 
              add = arg_checks)
  if(!inherits(gd[, "id"], "numeric") && !inherits(gd[, "id"], "character")){
    stop("The 'id' column in 'gd' must be of class numeric or character.")
  }
  assertIntegerish(gd[, "state"], add = arg_checks)
  assertNumeric(gd[, "time"], lower = 0, add = arg_checks)
  if(!missing(exact)){
    assertIntegerish(exact, lower = 1, upper = nrow(tmat), add = arg_checks)
  }
  assertIntegerish(maxit, lower = 0, add = arg_checks)
  assertNumeric(tol, add = arg_checks)
  assertLogical(verbose, add = arg_checks)
  
  # Check Spline arguments (See Eilers/Marx (2021) Section 2.3)
  
  #B-spline basis must be at least degree 0
  assertIntegerish(deg_splines, lower = 0, upper = Inf, max.len = 1, min.len = 1,
                   add = arg_checks)
  #Number of segments must be positive integer (deg_spline + 1 I think?)
  assertIntegerish(n_segments, lower = 1, upper = Inf, max.len = 1, min.len = 1,
                   add = arg_checks)
  #P-spline penalty order must be at least 1, must be smaller than
  #n_splines = deg_splines + n_segments
  assertIntegerish(ord_penalty, lower = 1, upper = deg_splines + n_segments -1,
                   min.len = 1, max.len = 1, add = arg_checks)
  
  assertFormula(formula, null.ok = TRUE, add = arg_checks)
  
  ## Data processing ---------------------------------------------------------

  #Retain original data (maybe scrap later)
  gd_orig <- gd
  
  #Change subject names to integers<<<<<<<<>>>>>>>>>>>>>>
  #Sort w.r.t. id and time
  gd <- gd[order(gd[, "id"], gd[, "time"]), ]
  original_subject_names <- unique(gd[, "id"]) #This can be numeric or character
  n_subjects <- length(original_subject_names)
  subject_names <- 1:n_subjects #New names, integer
  #Assign new names (now always integer)
  gd[, "id"] <- subject_names[match(gd[, "id"], original_subject_names)]
  
  #Transform data to matrix for faster operations
  gd <- as.matrix(gd)
  
  #Determine slices for each subject and store in hashmap (faster lookup)
  subject_slices <- new.env(parent = emptyenv(), size = n_subjects)
  for(subj in subject_names){
    subject_slices[[as.character(subj)]] <- which(gd[, "id"] == subj)
  }
  
  #Data transformations for probtrans compatibility
  #We shift all times by the minimum time, so that new minimum time becomes 0
  #because probtrans always calculates from 0 onward.
  min_time_orig <- min(gd[, "time"])
  gd[, "time"] <- gd[, "time"] - min_time_orig
  
  max_time <- max(gd[, "time"])
  
  tmat2 <- mstate::to.trans2(tmat)
  
  ## Derive parameters from data ---------------------------------------------

  n_states <- nrow(tmat) # no of states
  n_transitions <- nrow(tmat2) # no of transitions

  # Spline quantities
  n_splines <- deg_splines + n_segments
  
  #>>>>>>>>>>>>>>>>>>>>BIN LENGTH<<<<<<<<<<<<<<<<<<<<<
  #Determine length of bins and initialize them (using their mid and right-endpoints)
  #The bins are taken to be small enough that no two subsequent transitions 
  #of the same individual are present within any bin. 
  #To make it easy, we (for now) take the bin size to be the smallest 
  #difference between two consecutive observation times
  bin_length <- min(diff(sort(unique(gd[, "time"])))) #In article, represented by w
  
  bin_right <- seq(bin_length, max_time + bin_length, by = bin_length)
  n_bins <- length(bin_right) #Represented in article by U
  bin_middle <- (1:n_bins - 0.5)*bin_length #Represented by \overline{\tau_u}
  ## Derived Quantities Check ------------------------------------------------

  #If covariates specified.
  if(!missing(data)){
    subject_names_data <- unique(data[, "id"])
    
    #Check whether all observed subjects also have related covariates
    assertTRUE(all(subject_names %in% subject_names_data), add = arg_checks)
  }
  
  #Check if all states in gd are also in tmat
  assertTRUE(all(unique(gd[, "state"]) %in% 1:nrow(tmat2)), add = arg_checks)
  
  
  #Report assertions
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  
  
  ## Additional parameter determination ------------------------------------
  
  #Penalization matrix used in P-spline estimation
  penalization_matrix <- diff(diag(n_splines), diff=ord_penalty)

# EM Algorithm ------------------------------------------------------------

  
  #Pass by reference using environments (we don't want to copy objects all the time):
  #https://bookdown.org/content/d1e53ac9-28ce-472f-bc2c-f499f18264a3/reference.html
  
  #>>>>>>>>>>>>>>>Initiate initial estimates<<<<<<<<<<<<<<<<<<<
  #We need to initiate values for \theta = (\alpha, \beta)
  #with \alpha the B-spline coefficients 
  #and \beta the regression coefficients
  #We maximize over each transition separately,
  #so \vec{\alpha} can be represented using a n_splines x n_transitions matrix
  #with each column representing the estimates of \alpha_{gh} for that transition g->h
  
  EM_estimates <- new.env(parent = emptyenv())
  
  #>>>>>>>>>>>>>>>>>>>>>Initiate Spline coefficients<<<<<<<<<<<<<<<<<<<<<
  EM_estimates[["spline_coefficients_old"]] <- matrix(, nrow = n_splines, ncol = n_transitions)
  EM_estimates[["spline_coefficients_new"]] <- matrix(NA, nrow = n_splines, ncol = n_transitions)
  #>>>>>>>>>>>>>>>>>>>>>Initiate regression coefficients<<<<<<<<<<<<<<<<<<<<<
  
  #Initiate with infinite convergence criterion
  conv_criterion = Inf
  
  while(conv_criterion > tol){
    
    #We work separately for each transition
    for(transition in 1:n_transitions){
      ## E step ------------------------------------------------------------------
      
      #E.1
      #If data supplied: calculate risk-adjustment factors for each transition
      #for each subject: n_subjects x transitions matrix
      #Else: data <- NULL and in next functions adjust for this
      
      #E.2
      #Write probtrans_ODE or something, using ODE to:
      #Solve ODEs for all subjects (n), for all bins (U), determining:
      #P_{gh, i}(from, to). Maybe write only in 1 single direction? "forward"?
      
      
      
    
      # M Step ------------------------------------------------------------------
      
      #For this step, use JOPS::psPoisson() to obtain updated estimates
      #Remember that d_{gh,i}^u is a realisation of Poisson(Y_{g,i}^u exp(\eta_{gh,i}^u))
      #Check icpack for implementation of GLM Poisson regression. 
      #This will then yield coefficients of B-splines (for each transition)
      
      
      
    }
    ## Update estimates + convergence criterion ------------------
    
    
  }
    
  
  return(subject_slices)
  #DONT FORGET TO SHIFT times BACK BY 
  #min_time_orig
  
}