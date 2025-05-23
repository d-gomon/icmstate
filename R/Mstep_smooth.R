#'  Function for fitting proportional hazard model with baseline hazard
#'
#'
#'
#' @param fix_pars A list of fixed parameters during the EM procedure
#' @param EM_est A list of estimated quantities that change during the EM procedure
#' @param transno Transition number for the maximization procedure
#' @param from The state from which the relevant transition is taking place.
#' @param Pen Penalization matrix for the B-splines and regression coefficients.
#'
#'
#'
#'
#' @keywords internal
#' 
#' @return An object with fields:
#' \code{H} = hazards (matrix),
#' \code{cbx} = coefficient estimates (vector),
#' \code{lambda} = proposal for new lambda,
#' \code{ed} = effective dimension,
#' \code{G} = G matrix,
#' \code{ll} = log-likelihood,
#' \code{pen} = penalized part of log-likelihood,
#' \code{Mpen} = penalized M matrix
#'
#'
#'



Mstep_smooth <- function(fix_pars, EM_est, transno, from, Pen = Pen) {

  n_splines <- fix_pars[["n_splines"]] #=n_splines
  n_covariates <- fix_pars[["n_covariates"]] #=n_covariates
  use_RA <- fix_pars[["use_RA"]]
  B <- fix_pars[["Bspline_basis"]]
  X <- fix_pars[["mod_matrix"]]
  lambda <- EM_est[["lambda"]][transno]
  
  # Compute hazard, expected values and residuals
  coeff_splines <- c(EM_est[["coeff_old"]][1:n_splines, transno])
  if(!use_RA){
    coeff_covariates <- 0
  } else{
    coeff_covariates <- c(EM_est[["coeff_old"]][n_splines + (1:n_covariates), transno])
  }
  
  #browser()
  eta0 <- c(B %*% coeff_splines)  # log baseline hazard
  f <- c(X %*% coeff_covariates)  # covariate contribution to log hazard
  Eta <- outer(eta0, f, "+") #Basically add the two and stash them
  H <- exp(Eta) #Go to actual hazard from log-hazard
  Mu <- EM_est[["AtRisk"]][ , from, ] * H #Mu is equal to Y_{iu}*exp(\eta_{iu}), Basically hazard contribution times at risk prob
  Res <- EM_est[["NumTrans"]][ , transno, ] - Mu #Res is the y variable in article, the "Poisson outcomes"
  
  # Likelihood
  ll <- sum(EM_est[["NumTrans"]][ , transno, ] * Eta - Mu)
  #Likelihood that is being calculated here uses new EM_est values, 
  #but old coefficients for splines and covariates. This means 
  #the calculated likelihood will always be smaller than what we calculate in the E step
  pen <- lambda * t(coeff_splines) %*% Pen[1:n_splines, 1:n_splines] %*% coeff_splines/2
  
  # Construct equation system
  wb <- rowSums(Mu)
  wx <- colSums(Mu)
  qb <- t(B) %*% rowSums(Res)
  qx <- t(X) %*% colSums(Res)
  Mbb <- t(B) %*% diag(wb) %*% B         #diag.spam might not work appropriately
  Mxx <- t(X) %*% diag(wx) %*% X         #is diag() for sparse matrices.
  Mbx <- t(B) %*% Mu %*% X                    #Replace with diag() if issues arise
  M <- rbind(cbind(Mbb, Mbx), cbind(t(Mbx), Mxx))
  #browser()
  #Risk-adjustment case:
  if(use_RA){
    # Solve for new coefficients and check convergence
    Mpen <- M + Pen
    cnew <- solve(Mpen, c(qb, qx) + M %*% EM_est[["coeff_old"]][, transno])
    G <- solve(Mpen, M)
    ed <- sum(diag(G)[1:n_splines])
    EM_est[["coeff_new"]][, transno] <- cnew
  } else{ #No risk-adjustment, don't optimize the regression coefficient
    # Solve for new coefficients and check convergence
    Mpen <- M + Pen
    Mpen <- Mpen[1:n_splines, 1:n_splines]
    tempsol <- c(qb, qx) + M %*% EM_est[["coeff_old"]][, transno]
    tempsol <- tempsol[1:n_splines, , drop = FALSE]
    cnew <- solve(Mpen, tempsol)
    G <- solve(Mpen, M[1:n_splines, 1:n_splines])
    ed <- sum(diag(G)[1:n_splines])
    EM_est[["coeff_new"]][, transno] <- c(cnew,0)
  }
  
  if(!use_RA){ #If we don't have covariates, set "covariate" coefficient to 0
    EM_est[["coeff_new"]][n_splines + 1, transno] <- 0
  }
  u <- EM_est[["coeff_new"]][1:n_splines, transno]
  v <- t(u) %*% Pen[1:n_splines, 1:n_splines] %*% u / lambda
  # Because of error that Mpen was computationally singular
  EM_est[["lambda"]][transno] <- max(1e-4, min(1e+6, ed / v)) 
  
  #Update log-likelihood value
  EM_est[["complete_ll_new"]] <- ll
  #old output
  #output <- list(H = H, cbx = cbx, lambda = lambdanew, ed = ed, 
  #               G = G, ll = ll, pen = pen, Mpen = Mpen)
  return(EM_est)
}
