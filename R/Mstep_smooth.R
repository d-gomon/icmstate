#'  Function for fitting proportional hazard model with baseline hazard
#'
#' @import spam
#' @param Y Expected events (matrix)
#' @param R Expected risk sets (matrix)
#' @param X Covariates (matrix)
#' @param B B-spline basis
#' @param Pen Penalty matrix
#' @param lambda Smoothing parameter (number)
#' @param cbx Current coefficient estimates
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



Mstep_smooth <- function(fix_pars, EM_est, transno, from) {
  
  #function(Y, R, X, B, Pen, lambda, cbx)
  
  n_splines <- fix_pars[["n_splines"]] #=n_splines
  n_covariates <- fix_pars[["n_covariates"]] #=n_covariates
  B <- fix_pars[["Bspline_basis"]]
  lambda <- EM_est[["lambda"]][transno]
  
  # Compute hazard, expected values and residuals
  coeff_splines <- c(EM_est[["coeff_old"]][1:n_splines, transno])
  coeff_covariates <- c(EM_est[["coeff_old"]][n_splines + (1:n_covariates), transno])
  eta0 <- c(B %*% coeff_splines)  # log baseline hazard
  f <- c(X %*% coeff_covariates)  # covariate contribution to log hazard
  Eta <- outer(eta0, f, "+")
  H <- exp(Eta)
  Mu <- EM_est[["AtRisk"]][ , from, ] * H
  Res <- EM_est[["NumTrans"]][ , transno, ] - Mu
  
  # Likelihood
  ll <- sum(EM_est[["NumTrans"]][ , transno, ] * Eta - Mu)
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
  
  # Solve for new coeffcients and check convergence
  Mpen <- M + Pen
  cnew <- solve(Mpen, c(qb, qx) + M %*% EM_est[["coeff_old"]][, transno])
  G <- solve(Mpen, M)
  ed <- sum(diag(G)[1:n_splines])
  EM_est[["coeff_new"]][, transno] <- cnew
  u <- EM_est[["coeff_new"]][1:n_splines, transno]
  v <- t(u) %*% Pen[1:n_splines, 1:n_splines] %*% u / lambda
  # Because of error that Mpen was computationally singular
  EM_est[["lambda"]][transno] <- max(1e-4, min(1e+06, ed / v)) 
  
  #output <- list(H = H, cbx = cbx, lambda = lambdanew, ed = ed, 
  #               G = G, ll = ll, pen = pen, Mpen = Mpen)
}
