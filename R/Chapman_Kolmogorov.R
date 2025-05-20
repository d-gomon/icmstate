#' Chapman-Kolmogorov functions
#' 
#' @description
#' Functions which can be used to solver the Chapman-Kolmogorov equations. 
#' The functions are in the form as expected by deSolve:ode(). 
#' 
#' @returns Returns the derivative of the transition probability matrix P(s,t)
#'  with respect to time (forward: t, backward: s) as a list. 
#' 
#' @param t Time at which the derivative is required
#' @param state Values for the state in which the system resides at time t
#' (current "estimate" of transition matrix P). 
#' Must be a vector of length n_states * n_states containing the stacked 
#' columns of P: c(P_11, P_21, ... P_{n_states 1}, P_12, ..., P_{n_states n_states}).
#' @param parameters Parameters to derive the derivative. For P-splines, this is 
#' a list of coefficients, with each list entry (corresponding to a transition number)
#' containing a vector of n_splines coefficients.
#' @param n_states Number of states in the model.
#' 
#' @importFrom JOPS bbase
#' 
#' @keywords internal
#' 
#' 



# The Chapman-Kolmogorov functions
ChapKolm_fwd_smooth <- function(t, state, parms, fix_pars, subject) {
  #Extract some parameters
  n_states <- fix_pars[["n_states"]]
  n_covariates <- fix_pars[["n_covariates"]]
  n_splines <- fix_pars[["n_splines"]]
  max_time <- fix_pars[["max_time"]]
  n_segments <- fix_pars[["n_segments"]]
  deg_splines <- fix_pars[["deg_splines"]]
  tmat2_transids <- fix_pars[["tmat2_transids"]]
  use_RA <- fix_pars[["use_RA"]]
  #Current estimate of P: ode() will fill this in.
  P <- matrix(state, n_states, n_states)
  # Build dA matrix - containing intensities
  dA <- matrix(0, n_states, n_states)
  
  #pre-calculate some quantities
  coeffs_main <- parms[["coeff_old"]][1:n_splines, ]
  B <- bbase_singletime(t, xl = 0, xr = max_time, nseg = n_segments, 
                        bdeg = deg_splines)
  
  hazard_vec <- B %*% coeffs_main
  
  if(use_RA){
    ra_effect <- fix_pars[["mod_matrix"]][subject, ] %*% parms[["coeff_old"]][n_splines + (1:n_covariates), ]
    hazard_vec <- hazard_vec + ra_effect
  }
  
  hazard_vec <- exp(hazard_vec)
  dA[tmat2_transids] <- hazard_vec
  diag(dA) <- diag(dA) - rowSums(dA)
  # rate of change - Chapman-Kolmogorov equation
  dP <- P %*% dA
  # return the rate of change as list (requirement of deSolve:ode)
  list(c(dP)) #Returns the matrix, by columns
}


#' 
#' 
#' 
#' @keywords internal
#' @importFrom JOPS bbase
#' 

ChapKolm_bwd_smooth <- function(t, state, parms, fix_pars, subject) {
  #Extract some parameters  
  n_states <- fix_pars[["n_states"]]
  n_covariates <- fix_pars[["n_covariates"]]
  n_splines <- fix_pars[["n_splines"]]
  max_time <- fix_pars[["max_time"]]
  n_segments <- fix_pars[["n_segments"]]
  deg_splines <- fix_pars[["deg_splines"]]
  use_RA <- fix_pars[["use_RA"]]
  tmat2_transids <- fix_pars[["tmat2_transids"]]
  #Current estimate of P: ode() will fill this in.
  P <- matrix(state, n_states, n_states)
  # Build dA matrix
  dA <- matrix(0, n_states, n_states)
  
  #pre-calculate some quantities
  coeffs_main <- parms[["coeff_old"]][1:n_splines, ]
  B <- bbase_singletime(t, xl = 0, xr = max_time, nseg = n_segments, 
                        bdeg = deg_splines)
  
  hazard_vec <- B %*% coeffs_main
  
  if(use_RA){
    ra_effect <- fix_pars[["mod_matrix"]][subject, ] %*% parms[["coeff_old"]][n_splines + (1:n_covariates), ]
    hazard_vec <- hazard_vec + ra_effect
  }
  
  hazard_vec <- exp(hazard_vec)
  dA[tmat2_transids] <- hazard_vec
  
  diag(dA) <- diag(dA) - rowSums(dA)
  # rate of change - Chapman Kolmogorov equations
  dP <- dA %*% P
  # return the rate of change as list (requirement of deSolve:ode)
  list(c(dP))
}