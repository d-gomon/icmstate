#' Determine NPMLE for Multi State illness death Markov model using Frydman (1995)
#' 
#' @details For an illness death model (1 = healthy, 2 = ill, 3 = dead) estimate 
#' the NPMLE in the following form:
#'\describe{
#'   \item{\code{F12}:}{Cumulative distribution function of 1->2 transition;}
#'   \item{\code{F13}:}{Cumulative distribution function of 1->3 transition;}
#'   \item{\code{Lambda23}:}{Cumulative intensity of 2->3 transition;}
#' }
#' 
#' @param data A \code{data.frame} containing the columns named:
#' \describe{
#'   \item{\code{delta}:}{Did a transition from 1 -> 2 occur? (binary: 0 = no, 1 = yes); 
#'   In the left-truncated case, delta = 2 indicates initially observed in state 2.}
#'   \item{\code{Delta}:}{Was the transition to state 3 observed? (binary: 0 = no, 1 = yes);}
#'   \item{\code{L}:}{Left timepoint of interval censored transition to state 2 (numeric);}
#'   \item{\code{R}:}{Left timepoint of interval censored transition to state 2 (numeric);}
#'   \item{\code{time}:}{Time of event (transition to 3) or right-censoring in state 2 (numeric);}
#'   \item{\code{trunc}:}{(optional) Left-truncation time (numeric); Only used for entries with delta = 0 or 1.}
#' }
#' @param tol Tolerance of the EM algorithm. Algorithm will stop when the absolute difference 
#' between current mass estimates and new estimates is smaller than the tolerance
#' 
#' @references Frydman, H. (1995). Nonparametric Estimation of a Markov 
#' 'Illness-Death' Process from Interval- Censored Observations, with 
#' Application to Diabetes Survival Data. Biometrika, 82(4), 773-789. 
#' \doi{10.2307/2337344} 
#' 
#' @import checkmate
#' @export
#' 


msm_frydman <- function(data, tol = 1e-8){
  
  #Group 1: delta = 0, Delta = 0
  #Group 2: delta = 0, Delta = 1
  #Group 3: delta = 1, Delta = 0
  #Group 4: delta = 1, Delta = 1
  #Group 5: delta = 2, Delta = 0
  #Group 6: delta = 2, Delta = 1
  
  
  #-------------------Argument Checks-------------------------------#
  arg_checks <- makeAssertCollection()
  assertDataFrame(data, min.cols = 5, max.cols = 6, add = arg_checks)
  data_ncols <- ncol(data)
  if(data_ncols == 5){
    assertNames(names(data), must.include = c("delta", "Delta", "L", "R", "time"), 
                add = arg_checks)  
  } else if(data_ncols == 6){
    assertNames(names(data), must.include = c("delta", "Delta", "L", "R", "time", "trunc"), 
                add = arg_checks)
    assertNumeric(data[["trunc"]], lower = 0, any.missing = FALSE, add = arg_checks)
    #TODO: If missing truncation time, set to NA
  }
  assertSubset(data[["delta"]], c(0,1,2), add = arg_checks)
  assertSubset(data[["Delta"]], c(0,1), add = arg_checks)
  assertNumeric(data[["L"]], lower = 0, add = arg_checks)
  assertNumeric(data[["R"]], lower = 0, add = arg_checks)
  assertNumeric(data[["time"]], lower = 0, any.missing = FALSE, add = arg_checks)
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)

  
  #-----------------Data transformations------------------------#
  ltrunc <- FALSE
  #Are we working in the left-truncated case?
  if(data_ncols == 6){
    ltrunc <- TRUE
  }
  #Re-arrange so that columns are always in the following order:
  col_order <- c("delta", "Delta", "L", "R", "time")
  if(isTRUE(ltrunc)){
    col_order <- c(col_order, "trunc")
  }
  data <- data[, col_order]
  #Transform to matrix for faster operations
  matdata = as.matrix(data)
  #Retain row names
  rownames(matdata) <- rownames(data)
  #Determine indices of groups:
  Group1_idx = which(matdata[, 1] == 0 & matdata[, 2] == 0)
  Group2_idx = which(matdata[, 1] == 0 & matdata[, 2] == 1)
  Group3_idx = which(matdata[, 1] == 1 & matdata[, 2] == 0)
  Group4_idx = which(matdata[, 1] == 1 & matdata[, 2] == 1)
  Group5_idx = which(matdata[, 1] == 2 & matdata[, 2] == 0)
  Group6_idx = which(matdata[, 1] == 2 & matdata[, 2] == 1)
  
  #Assertions when delta == 1: L and R both non missing and L < R, time > R
  #Need to check whether L and R non-missing in Groups 3 and 4 and whether
  #L < R (strictly) otherwise we can get errors further on in code
  assert(!anyMissing(matdata[c(Group3_idx, Group4_idx), c(3,4)]))
  assert(all(matdata[c(Group3_idx, Group4_idx), 3] < matdata[c(Group3_idx, Group4_idx), 4]), 
         .var.name = "Left bound L in interval censoring must be smaller than right bound R.")
  #Assertion for time > R (Groups 4 & 6)
  assert(all(matdata[c(Group4_idx, Group6_idx), 4] <= matdata[c(Group4_idx, Group6_idx), 5]),
  .var.name = "Right bound R must be smaller than 
         event time.")
  if(isTRUE(ltrunc)){
    #Assertions for truncation:
    assert(all(matdata[c(Group1_idx, Group2_idx), 6] <= matdata[c(Group1_idx, Group2_idx), 5]), 
           .var.name = "Truncation time must be smaller/equal than event time.")
    assert(all(matdata[c(Group3_idx, Group4_idx), 6] <= matdata[c(Group3_idx, Group4_idx), 3]), 
           .var.name = "Truncation time must be smaller/equal than left bound L in interval censoring.")
    assert(all(matdata[c(Group5_idx, Group6_idx), 6] <= matdata[c(Group5_idx, Group6_idx), 5]), 
           .var.name = "Truncation time must be smaller/equal than right bound R in interval censoring.")
  }
  

  #Set interval censoring times to missing when transition hasn't happened
  matdata[c(Group1_idx, Group2_idx), c(3, 4)] <- NA
  #Set left bound interval censoring to missing when unknown (Group 5, 6)
  matdata[c(Group5_idx, Group6_idx), 3] <- NA
  data_idx <- list(
    #Data + group data
    matdata = matdata,
    Group1_idx = Group1_idx,
    Group2_idx = Group2_idx,
    Group3_idx = Group3_idx,
    Group4_idx = Group4_idx,
    #Length indicators for the groups
    J = length(Group1_idx),
    Ktilde = length(Group2_idx),
    Ntilde = length(Group4_idx),
    gr3tilde = length(Group3_idx),
    M = length(Group3_idx) + length(Group4_idx),
    N_star = nrow(matdata),
    #Subset groups to separate data frames
    matdata_g1 = matdata[Group1_idx, , drop = FALSE],
    matdata_g2 = matdata[Group2_idx, , drop = FALSE],
    matdata_g3 = matdata[Group3_idx, , drop = FALSE],
    matdata_g4 = matdata[Group4_idx, , drop = FALSE],
    matdata_g34 = matdata[c(Group3_idx, Group4_idx), , drop = FALSE],
    ltrunc = ltrunc
  )
  
  

  #Failure times and their lengths
  t_n_star = sort(unique(data_idx$matdata_g4[, 5]))
  d_n = table(data_idx$matdata_g4[, 5])
  N = length(t_n_star)
  e_k_star <- sort(unique(data_idx$matdata_g2[, 5]))
  c_k = table(data_idx$matdata_g2[, 5])
  K = length(e_k_star)
  
  data_idx <- c(data_idx, list(
    t_n_star = t_n_star,
    d_n = d_n,
    N = N,
    e_k_star = e_k_star,
    c_k = c_k,
    K = K
  ))
  
  
  #Some truncation data, add to data_idx
  if(isTRUE(ltrunc)){
    matdata_g5 = matdata[Group5_idx, , drop = FALSE]
    matdata_g6 = matdata[Group6_idx, , drop = FALSE]
    t_n_star_onlytrunc = sort(unique(matdata_g6[, 5]))
    d_n_onlytrunc <- table(matdata_g6[, 5])
    B <- length(t_n_star_onlytrunc)
    Btilde <- length(Group6_idx)
    t_n_star_trunc <- sort(unique(c(t_n_star, t_n_star_onlytrunc)))
    d_n_star_trunc <- table(c(data_idx$matdata_g4[, 5], matdata_g6[, 5]))
    v_r <- matdata[c(Group1_idx, Group2_idx, Group3_idx, Group4_idx), 6]
    
    data_idx <- c(data_idx, list(
      matdata_g5 = matdata_g5,
      matdata_g6 = matdata_g6,
      t_n_star_onlytrunc = t_n_star_onlytrunc,
      d_n_onlytrunc = d_n_onlytrunc,
      B = B,
      Btilde = Btilde,
      t_n_star_trunc = t_n_star_trunc,
      d_n_star_trunc = d_n_star_trunc,
      v_r = v_r,
      Z = length(Group5_idx) + length(Group6_idx)
    ))
  }

  
  #----------------Find support of CDF--------------------------#
  #See Frydman (1995) Page 777#
  #Output: sets {L_tilde}, {R_tilde} and {Q_i}
  supportMSM <- support_frydman(data_idx)  
  #Returns list with A, unique_L_times, L_tilde, R_tilde, Q_mat, I.

  
  #---------------EM algorithm---------------------------------#
  #Initiate values for z and lambda (take trivial values for now)
  z_init <- rep(1, supportMSM$I + K)/(supportMSM$I+K)
  lambda_init <- rep(0.5, data_idx$N)
  
  #output: values of z and lambda which maximize likelihood
  if(isTRUE(ltrunc)){
    z_lambda <- EM_solver_trunc(data_idx = data_idx, supportMSM = supportMSM,
                          z = z_init, lambda = lambda_init, tol = tol) 
  } else{
    z_lambda <- EM_solver(data_idx = data_idx, supportMSM = supportMSM,
                          z = z_init, lambda = lambda_init, tol = tol)  
  }

  
  
  #Use estimated values of z and lambda to calculate cdf of transitions.
  cdf_survival <- cdf_calculator(z_lambda = z_lambda, supportMSM = supportMSM,
                                 data_idx = data_idx)
  out <- list(data_idx = data_idx, supportMSM = supportMSM,
              z_lambda = z_lambda, cdf = cdf_survival)
  class(out) <- "msmFrydman"
  return(out)
}



















