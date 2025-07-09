#' Smooth hazard estimation for general multi-state model with interval censored
#' transitions 
#' 
#' @description For a general Markov chain multi-state model with interval censored 
#' transitions smoothly estimate the (log-)hazard (see Details). The estimates 
#' are also returned as an \code{\link[mstate:msfit]{msfit}} object. The smallest time 
#' in the data will be set to zero.
#' Note that the hazard is modelled on LOG-scale!! This means that the B-splines 
#' will represent the log-hazard. For underlying weibull hazard, this means that
#' the hazard will always be a first degree polynomial in log(time). Keep this in 
#' mind when choosing the degree of the B-splines.
#' 
#' 
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param tmat A transition matrix as created by \code{\link[mstate:transMat]{transMat}}
#' @param formula Formula to interpret in data for covariates.
#' @param data A \code{data.frame} containing a column called \code{'id'} 
#' (identifying the subjects in \code{gd}) and 
#' variables which to interpret in \code{formula}.
#' @param deg_splines Degree to use for the B-spline basis functions. Defaults 
#' to 3 (cubic B-splines).
#' @param n_segments Number of segments to use for the P-splines. The
#' segments will space the domain evenly. According to Eilers & Marx (2021), it 
#' it OK to choose this number very large. Default = 20. 
#' @param ord_penalty Order of the P-spline penalty (penalty on the difference 
#' between d-order differences of spline coefficients). See Eilers & Marx 
#' Section 2.3. Defaults to 2.
#' @param n_bins Number of bins to use for estimation. The higher this number, the 
#' higher the 'resolution' of the estimates (see Details). Guideline is to have the bin size 
#' smaller than the smallest 2 consecutive observation times for any subject.
#' Large values increase computational load. Default = 100.
#' @param maxit Maximum number of iterations. Default = 100.
#' @param tol Tolerance of the convergence procedure in the E-step. A change in the value of 
#' the observed-data likelihood in an iteration of less than \code{tol} will make the procedure stop.
#' @param Mtol Tolerance of the convergence procedure of the M-step. The M-step 
#' consists of an Iteratively Reweighted Least Squares (IRLS) procedure, where 
#' the (unobserved) complete-data likelihood is maximized. Default is \code{1e-4}.
#' @param verbose Should iteration messages be printed? Default is FALSE
#' @param ode_solver The integrator to use for solving the ODE's. See 
#' \code{\link[deSolve:ode]{ode()}}. By default, the "lsoda" solver will be used.
#' @param ridge_penalty The ridge penalty to use for estimating risk-adjustment 
#' coefficients. Default = 1e-06.
#' 
#' 
#' 
#' @details
#' This function estimates (for each transition in the MSM) the log hazard
#' at the midpoint of bins \eqn{\eta_{iu}}{eta[iu]} given by:
#' \deqn{\log(h_{iu}) = \eta_{iu} = \sum_{m=1}^{M} b_{um} \alpha_{m} + \sum_{j=1}^{p} x_{ij} \beta_{j}}{eta[iu] = sum(b[um] * alpha[m], m=1..M) + sum(x[ij] * beta[j], j=1..p)}
#' with \eqn{M} the number of splines (given by \code{n_segments + deg_splines}). u 
#' represents a single bin \eqn{(\tau_{u-1}, \tau_u]} (total number specified by \code{n_bins}) such that 
#' \eqn{b_{um} = B_m(\overline{\tau}_u)} with \eqn{\overline{\tau}_u} the midpoint 
#' of bin u. Each spline gets assigned a single parameter \eqn{\alpha_m}, which 
#' can heuristically be thought of as the approximate value of the spline. 
#' All bins have the same length and cover the total observed 
#' time-period evenly. The second part of the equation corresponds to the 
#' subject specific covariates and coefficients respectively. The function returns 
#' values for \eqn{\theta = (\alpha, \beta)}, thereby allowing users to reconstruct 
#' the (log-)hazard.
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
#' @importFrom utils hashtab sethash
#' @importFrom stats model.matrix
#' @import checkmate
#' @export
#' 
#' 
#' @examples
#' #Generate data from illness-death model with Weibull hazards.
#' wshapes <- c(0.5,0.5,2)
#' wscales <- c(5, 10, 10/gamma(1.5))
#' ID_Weib <- sim_weibmsm(tmat = trans.illdeath(), shape = wshapes,
#'                        scale = wscales, n_subj = 80, obs_pars = c(2, 0.5, 20), 
#'                        startprobs = c(0.9, 0.1, 0))
#' 
#' smoothID_Weib <- smoothmsm(gd = ID_Weib, tmat = trans.illdeath(), ord_penalty = 1,
#' deg_splines = 1, tol = 1e-5, maxit = 20)
#' 
#' #Let's plot the estimates against the truth
#' x <- seq(0, max(ID_Weib$time), 0.01)
#' #True hazard (see pweibull details)
#' haz1 <- -pweibull(x, wshapes[1], wscales[1], lower = FALSE, log = TRUE)
#' haz2 <- -pweibull(x, wshapes[2], wscales[2], lower = FALSE, log = TRUE)
#' haz3 <- -pweibull(x, wshapes[3], wscales[3], lower = FALSE, log = TRUE)
#' 
#' plot(smoothID_Weib$smoothmsfit, main = "SMOOTH (dashed = true)")
#' lines(x, haz1, col = "black", lwd = 2, lty = 2)
#' lines(x, haz2, col = "red", lwd = 2, lty = 2)
#' lines(x, haz3, col = "green", lwd = 2, lty = 2)
#' 
#' 
#' 



smoothmsm <- function(gd, tmat, formula, data,
                      deg_splines = 3, n_segments = 20, ord_penalty = 2, n_bins = 100,
                      maxit = 100, tol = 1e-4, Mtol = 1e-4, 
                      verbose = FALSE, ode_solver = "lsoda",
                      ridge_penalty = 1e-06){
  
# Pre-processing ---------------------------------------

  if(missing(formula)){
    formula <- NULL
  }

  ## Argument checks ---------------------------------------------------------
  call <- match.call()
  #conv_crit <- match.arg(conv_crit, choices = c("haz", "prob", "lik"))
  
  arg_checks <- makeAssertCollection()
  assertMatrix(tmat, add = arg_checks)
  assertDataFrame(gd, min.cols = 3, max.cols = 3, add = arg_checks)
  assertNames(names(gd), must.include = c("id", "state", "time"), 
              add = arg_checks)
  if(!inherits(gd[, "id"], "numeric") && !inherits(gd[, "id"], "character") 
     && !inherits(gd[, "id"], "integer")){
    stop("The 'id' column in 'gd' must be of class numeric or character.")
  }
  assertIntegerish(gd[, "state"], add = arg_checks)
  assertNumeric(gd[, "time"], lower = 0, add = arg_checks)
  #if(!missing(exact)){
  #  assertIntegerish(exact, lower = 1, upper = nrow(tmat), add = arg_checks)
  #}
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

  #Retain original data on states (maybe scrap later)
  gd_orig <- gd
  
  #Change subject names to integers<<<<<<<<>>>>>>>>>>>>>>
  #Sort w.r.t. id and time
  gd <- gd[order(gd[, "id"], gd[, "time"]), ]
  original_subject_names <- unique(gd[, "id"]) #This can be numeric or character
  n_subjects <- length(original_subject_names)
  
  subject_names <- 1:n_subjects #New names, integer
  #If covariates specified.
  if(!missing(data)){
    data <- data[order(data[, "id"]), ]
    if(nrow(data) != n_subjects){
      stop("Number of subjects in 'gd' and 'data' is not equal.")
    }
    subject_names_data <- unique(data[, "id"])
    #Check whether all observed subjects also have related covariates
    assertTRUE(all(original_subject_names %in% subject_names_data))
    #Extract only relevant entries
    #Translate data names to new names
    data[, "id"] <- subject_names[match(data[, "id"], original_subject_names)]
  }
  gd[, "id"] <- subject_names[match(gd[, "id"], original_subject_names)]
  
  #Transform data to matrix for faster operations
  gd <- as.matrix(gd)
  
  
  #Data transformations for probtrans compatibility
  #We shift all times by the minimum time, so that new minimum time becomes 0
  #because probtrans always calculates from 0 onward.
  min_time_orig <- min(gd[, "time"])
  gd[, "time"] <- gd[, "time"] - min_time_orig
  
  
  
  tmat2 <- mstate::to.trans2(tmat)
  tmat2_transids <- as.matrix(tmat2[, c("from", "to")])
  ## Derive parameters from data ---------------------------------------------

  n_states <- nrow(tmat) # no of states
  n_transitions <- nrow(tmat2) # no of transitions
  
  # Spline quantities
  n_splines <- deg_splines + n_segments
  
  #>>>>>>>>>>>>>>>>>>>>BIN LENGTH<<<<<<<<<<<<<<<<<<<<<
  
  #Switched to specifying n_bins in advance, requiring the calculation of bin_length
  #from the OriGinal time frame. We extend the time-frame a little to compensate.
  max_time_og <- max(gd[, "time"]) * 1.05
  bin_length <- max_time_og/n_bins
  #Rescale time to integers on bin scale (so (0, 1] is bin 1, (1, 2] is bin 2 etc...)
  gd[, "time"] <- gd[, "time"]/bin_length
  max_time <- n_bins
  bin_right <- 1:max_time
  bin_middle <- (1:n_bins) - 0.5
  
  
  #Old method to determine bin length and n_bins
  #Determine length of bins and initialize them (using their mid and right-endpoints)
  #The bins are taken to be small enough that no two subsequent transitions 
  #of the same individual are present within any bin. 
  #To make it easy, we (for now) take the bin size to be the smallest 
  #difference between two consecutive observation times of single subject.
  #diff_times <- diff(unique(gd[, "time"]))
  #bin_length <- min(diff_times[diff_times > 0]) #In article, represented by w
  
  #Scale data by bin length (so bin endpoints are always at 1,2,...)
  #gd[, "time"] <- gd[, "time"]/bin_length
  #max_time <- ceiling(max(gd[, "time"])) + 1L
  #bin_right <- 1:max_time
  #n_bins <- max_time #Represented in article by U
  #bin_middle <- (1:n_bins) - 0.5 #Represented by \overline{\tau_u}
  
  ## Derived Quantities Check ------------------------------------------------
  
  #Check if all states in gd are also in tmat
  assertTRUE(all(unique(gd[, "state"]) %in% 1:nrow(tmat)), add = arg_checks)
  
  #Report assertions
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)
  
  
  ## Additional parameter determination ------------------------------------
  
  #Penalization matrix used in P-spline estimation
  penalization_matrix <- diff(diag(n_splines), diff = ord_penalty)
  
  
  #B-spline basis matrix used in M-step
  Bspline_basis <- JOPS::bbase(bin_middle, xl = 0, xr = max_time, 
                               nseg = n_segments, bdeg = deg_splines)
  
  ## Risk adjustment ----------------------------------------------------
  if(!missing(data) & !missing(formula)){
    rhs <- formula[-2]
    mod_matrix <- model.matrix(object = rhs, data = data)
    #Remove intercept column
    intercept_pos <- match("(Intercept)", colnames(mod_matrix))
    if(!is.na(intercept_pos)){
      mod_matrix <- mod_matrix[, -intercept_pos, drop = FALSE]
    }
    n_covariates = ncol(mod_matrix)
    use_RA <- TRUE
  } else{
    #Make an zero-dummy model matrix
    mod_matrix <- matrix(0, nrow = n_subjects, ncol = 1)
    colnames(mod_matrix) <- "Dummy"
    n_covariates = 1
    use_RA <- FALSE
  }
  
  
  # Ridge Regularization
  if(use_RA){
    ridge_penalty <- ridge_penalty #User could choose this...  
  } else{
    ridge_penalty <- 0 #No ridge penalty when no risk-adjustment
  }
  
  D <- diff(diag(n_splines), diff = ord_penalty)
  Pdiff = Pridge = 0 * diag(n_splines + n_covariates)
  Pdiff[1:n_splines, 1:n_splines] = t(D) %*% D   #I don't understand why, but sure.
  diag(Pridge)[n_splines + (1:n_covariates)] = ridge_penalty
  
  ## bbase - caching ---------------------------------------------------------
  
  #We can cache some quantities which we need to calculate bases.
  #Majorly improved computation speed!
  bbase_dx <- (max_time - 0)/n_segments
  bbase_n_knots <- 2*deg_splines + 2
  bbase_t_D <- diff_D_t(diag(bbase_n_knots), differences = deg_splines + 1) / 
    (gamma(deg_splines + 1) * bbase_dx ^ deg_splines)
  bbase_t_D_trunc <- diff_D_t(diag(bbase_n_knots - 1), differences = deg_splines + 1) / 
    (gamma(deg_splines + 1) * bbase_dx ^ deg_splines)
  bbase_const <- (-1)^(deg_splines + 1)
  
  
  ## Fixed parameters - single list -----------------------
  #We could make this into an environment if we run into performance issues
  #Using hashtab, simply perform sethash() multiple times, with character names "n_segments" etc.
  fix_pars <- list(gd = gd,
                   deg_splines = deg_splines,
                   n_segments = n_segments,
                   ord_penalty = ord_penalty,
                   n_bins = n_bins,
                   maxit = maxit,
                   tol = tol,
                   #prob_tol = prob_tol,
                   n_states = n_states,
                   n_transitions = n_transitions,
                   n_splines = n_splines,
                   bin_length = bin_length,
                   penalization_matrix = penalization_matrix,
                   n_covariates = n_covariates,
                   n_coefficients = n_splines + n_covariates,
                   ode_solver = ode_solver,
                   n_subjects = n_subjects,
                   max_time = max_time,
                   tmat = tmat,
                   tmat2 = tmat2,
                   tmat2_transids = tmat2_transids,
                   Bspline_basis = Bspline_basis,
                   use_RA = use_RA,
                   mod_matrix = mod_matrix,
                   bbase_dx = bbase_dx,
                   bbase_t_D = bbase_t_D,
                   bbase_t_D_trunc = bbase_t_D_trunc,
                   bbase_const = bbase_const)
  
  
  #Subjects are now always named 1:n_subjects
  #Original subject names can be recovered from original_subject_names
  
  #max_time is the maximum time in the scaled timeframe, representing 
  #the right endpoint of the final "bin"
  
  
  #Determine slices for each subject and store in hashmap (faster lookup)
  subject_slices <- hashtab(type = "identical", size = n_subjects)
  for(subj in subject_names){
    sethash(subject_slices, subj, which(gd[, "id"] == subj))
  }
  #Takes a long time initially, but compensates a lot later.
  
  
  ## Check if n_bins chosen appropriately -----------------------
  min_diff_time_subj <- Inf
  for(subj in subject_names){
    #Smallest time difference between 3 consecutive observations
    temp_min_diff <- min(diff(qwe, lag = 2))
    if(temp_min_diff < min_diff_time_subj){
      min_diff_time_subj <- temp_min_diff
    }
  }
  
  if(bin_length > min_diff_time_subj){
    n_bins_suggested <- max_time_og/min_diff_time_subj
    warning(paste0("Chosen 'n_bins' does not appropriately cover the data. We suggest at least ", max(ceiling(n_bins_suggested), n_splines), " n_bins."))
  }

  
# EM Algorithm ------------------------------------------------------------

  
  #Pass by reference using environments (we don't want to copy objects all the time):
  #https://bookdown.org/content/d1e53ac9-28ce-472f-bc2c-f499f18264a3/reference.html
  
  
  ## EM estimates list ----------------------------------------------
  #We need to initiate values for \theta = (\alpha, \beta)
  #with \alpha the B-spline coefficients 
  #and \beta the regression coefficients
  #We maximize over each transition separately,
  #so c(\vec{\alpha}, \vec{\beta}) can be represented using a 
  #(n_splines + n_coefficients) x n_transitions matrix
  #with each column representing the estimates of the coefficients for that transition g->h
  
  #EM_est <- hashtab()
  #hashtab() is slower than just using list when accessing often :(
  EM_est <- list()
  
  #>>>>>>>>>>>>>>>>>>>>>Initiate Spline + regression coefficients<<<<<<<<<<<<<<<<<<<<<
  
  # We stash the coefficients as follows \theta = c(\alpha, \beta)
  EM_est[["coeff_old"]] <- matrix(c(rep(log(1/n_bins), n_splines), 
                                           rep(0, n_covariates)), 
                                         nrow = n_splines + n_covariates, 
                                         ncol = n_transitions)
  EM_est[["coeff_new"]] <- matrix(c(rep(log(1/n_bins), n_splines), 
                                           rep(0, n_covariates)), 
                                         nrow = n_splines + n_covariates, 
                                         ncol = n_transitions)
  
  #Maybe think of other initial conditions
  
  #>>>>>>>>>>>>>>>>>>>>>>Initiate expected value matrices<<<<<<<<<<<<<<<<<<<
  EM_est[["AtRisk"]] <- array(0, dim = c(n_bins, n_states, 
                                                       n_subjects))
  EM_est[["NumTrans"]] <- array(0, dim = c(n_bins, n_transitions, 
                                                       n_subjects))
  
  
  #>>>>>>>>>>>>>>>>>>>>>>>Keep track of likelihood<<<<<<<<<<<<<<<<<<<<<<<<<<
  
  # Observed log likelihood (calculated in E-step)
  EM_est[["loglik_old"]] <- -Inf
  EM_est[["loglik_new"]] <- -Inf
  
  #"Complete" (unobserved) Poisson log likelihood (calculated in M-step)
  EM_est[["complete_ll_old"]] <- rep(-Inf, n_transitions)
  EM_est[["complete_ll_new"]] <- rep(-Inf, n_transitions)
  
  
  #>>>>>>>>>>>>>>>>>>>>>>Penalization coefficient lambda<<<<<<<<<<<<<<<<<<<<
  # Separate penalization coefficient for each transition.
  EM_est[["lambda"]] <- rep(100, n_transitions)
  
  #Initiate with infinite convergence criterion
  conv_criterion = Inf
  
  ll_history <- vector(mode = "numeric", length = maxit)
  it_num <- 1
  while(conv_criterion > tol & it_num < maxit + 1){
    cat(it_num)
    ## E step ------------------------------------------------------------------
    
    #cat("Start Estep iteration", it_num, "\n")
    #E.2
    #Write probtrans_ODE or something, using ODE to:
    #Solve ODEs for all subjects (n), for all bins (U), determining:
    #P_{gh, i}(from, to).
    
    #Called for "side" effects: Updates the values in EM_est (hashtable)
    #In particular: "AtRisk" and "NumTrans"
    EM_est <- Estep_smooth(fix_pars = fix_pars, subject_slices = subject_slices,
                 EM_est = EM_est, it_num = it_num)
    
    #cat("E-ll using estimates from it ", it_num-2, ": ", EM_est[["loglik_old"]], "\n")
    #cat("E-ll using estimates from it ", it_num-1, ": ", EM_est[["loglik_new"]], "\n")
    
    
    #We don't want to update lambda in every iteration, only when we have converged
    #lambda_init <- EM_est[["lambda"]]
    Mstep_conv_criterion = Inf
    Mtol <- Mtol
    #We work separately for each transition in the M step
    while(Mstep_conv_criterion > Mtol){ #How many Mstep iterations do we perform?
      Mstep_ll <- 0
      #EM_est[["lambda"]] <- lambda_init #If we do another iteration, we should use original lambda value
      for(transno in 1:n_transitions){
        #Update coefficients if we want to do another iteration of M-step
        EM_est[["coeff_old"]] <- EM_est[["coeff_new"]]
        ## M Step ------------------------------------------------------------------
        #For this step, the idea is to use a variant of JOPS::psPoisson() to obtain updated estimates
        #Remember that d_{gh,i}^u is a realisation of Poisson(Y_{g,i}^u exp(\eta_{gh,i}^u))
        #This will then yield coefficients of B-splines (for each transition)
        #and due to the manual adjustment of psPoisson() in Mstep_smooth
        #also the regression coefficients
        
        Pen <- EM_est[["lambda"]][transno] * Pdiff + Pridge
        from <- tmat2[transno, "from"]
        to <- tmat2[transno, "to"]
        #Called for side-effects, updates EM_est: "lambda" and "coeff_old" and "coeff_new"
        EM_est <- Mstep_smooth(fix_pars = fix_pars, EM_est = EM_est, transno = transno, from = from, Pen = Pen)
        #returns log-likelihood contribution for transno, using old coefficient estimates
        
        
        #cat("Start Mstep transition", transno, " completed. Log-likelihood contribution ", Mstep_ll_temp, "\n")
      }
      
      Mstep_conv_criterion <- sum(EM_est[["complete_ll_new"]]) - sum(EM_est[["complete_ll_old"]])
      EM_est[["complete_ll_old"]] <- EM_est[["complete_ll_new"]]
      if(Mstep_conv_criterion < 0){ #Sometimes M step can go in the negative direction
        Mstep_conv_criterion = Inf
      }
    }
    #cat("Mstep likelihood: ", sum(EM_est[["complete_ll_new"]]), "Mstep diff: ", 
    #    Mstep_conv_criterion, "\n")
    
    ## Check convergence  ------------------
    ll_history[it_num] <- EM_est[["loglik_old"]]
    ll_dif <- EM_est[["loglik_new"]] - EM_est[["loglik_old"]]
    conv_criterion <- ll_dif
    if(conv_criterion < 0){ #Sometimes we can have a decrease in likelihood...
      conv_criterion = Inf
    }
    
    cat("End of Iteration: ",it_num, "E-ll: ", EM_est[["loglik_new"]], 
        "ll_diff: ", ll_dif, "\n")
    
    
    ## Update EM estimates -------------------
    
    EM_est[["coeff_old"]] <- EM_est[["coeff_new"]]
    EM_est[["loglik_old"]] <- EM_est[["loglik_new"]]
    
    it_num <- it_num + 1
  }

# Output Creation ---------------------------------------------------------

  #Create a smoothmsfit object which we will use for plotting
  #Here we also shift back time to original frame:
  #First, scale time back by * bin_length
  #Then shift back time by: +min_time_orig
  smoothmsfit <- create_smoothmsfit(fix_pars = fix_pars, EM_est = EM_est)
  coeff_out <- EM_est[["coeff_old"]]
  rownames(coeff_out) <- c(paste0("B", 1:n_splines), colnames(fix_pars[["mod_matrix"]]))
  colnames(coeff_out) <- paste0("trans", 1:n_transitions)
  out <- list(log_coeff = coeff_out,
              AtRisk = EM_est[["AtRisk"]],
              NumTrans = EM_est[["NumTrans"]],
              loglik = ll_history,
              fix_pars = fix_pars,
              EM_est = EM_est,
              min_time = min_time_orig,
              smoothmsfit = smoothmsfit)

  #Assign 'smoothmsm' class so we can use 'smoothmsm' functions on output
  class(out) <- "smoothmsm"
  return(out)
}



#' Create 'smoothmsfit' object from 'smoothmsm' output
#' 
#' The goal is to create a data.frame which contains all the necessary 
#' information for plotting hazards, allowing for risk-adjusted curves as well.
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' 
#' 
#' 

create_smoothmsfit <- function(fix_pars, EM_est){
  #We create a 'msfit' part of this function containing the usual
  #$Haz ($time, $Haz, $trans), $varHaz, $trans components for the baseline hazard
  #Also, with Bsplines we actually recover the (non-cumulative) hazard, which we store as
  #$haz ($time, $haz, $trans), $varhaz, $trans
  #Additionally, we will have $mod_matrix (containing the covariates for subjects)
  #and $BsplineBasis for the basis used 
  
  out <- vector(mode = "list")
  out$mod_matrix <- fix_pars[["mod_matrix"]]
  out$Bspline_basis <- fix_pars[["Bspline_basis"]]
  #Coefficients will be a list containing the $spline and $covariate coefficients
  #Each of these will be a matrix of size n_splines or n_covariates times n_transitions
  out$coefficients <- list(spline = EM_est[["coeff_old"]][1:fix_pars[["n_splines"]], ],
                           covariate = EM_est[["coeff_old"]][fix_pars[["n_splines"]] + (1:fix_pars[["n_covariates"]]), ])
  out$trans <- fix_pars[["tmat"]]
  #Now we add the 'msfit' part (baseline hazard)
  #For now, we can't calculate the variances yet.
  #haz is a n_times x n_transitions matrix, containing log hazards at each time for each transition
  haz <- fix_pars[["Bspline_basis"]] %*% out[["coefficients"]][["spline"]]
  #Obtain hazard function
  out$haz <- data.frame(time = rep((1:fix_pars[["max_time"]]) * fix_pars[["bin_length"]], 
                             fix_pars[["n_transitions"]]),
                  haz = exp(as.vector(haz)),
                  trans = rep(1:fix_pars[["n_transitions"]], each = fix_pars[["max_time"]]))
  #Obtain cumulative hazard function (can be fit using plot.msfit)
  out$Haz <- data.frame(time = out[["haz"]][["time"]],
                  Haz = as.vector(apply(exp(haz), 2, cumsum)),
                  trans = out[["haz"]][["trans"]])
  
  #Remove "msfit" later. Make plot.smoothmsfit which allows to just use plot.msfit
  class(out) <- c("msfit", "smoothmsfit")
  return(out)
}

