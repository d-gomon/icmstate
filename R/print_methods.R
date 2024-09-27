#' @title Print a "npmsm" object
#' 
#' @description Print some details of a \code{\link{npmsm}} fit
#' 
#' @param x An object of class \code{"npmsm"}
#' @param ... Additional arguments to print
#' 
#' @return A summary of the fitted model will be displayed in the console
#'
#' @export
#'



print.npmsm <- function(x, ...){
  method_used <- ifelse(x$method == "multinomial", "Multinomial", "Poisson")
  crit_used <- switch(x$conv_crit, 
                      "haz" = "change in estimated intensities",
                      "prob" = "change in estimated probabilities",
                      "lik" = "change in observed-data likelihood")
  cat(paste0(method_used, " complete-data likelihood maximized through EM algorithm.\n"))
  cat(paste0("Convergence criterion: ", crit_used, " < ", x$tol, "\n"))
  if(x$it < x$maxit){
    cat(paste0("Algorithm converged after ", x$it, " out of max ", x$maxit, " iterations.\n"))  
  } else{
    cat(paste0("Algorithm failed to converge after ", min(x$it, x$maxit), " out of max ", x$maxit, " iterations.\n"))
    cat("Consider increasing 'maxit'.\n")
  }
  cat(paste0("Log likelihood value: ", x$ll, "\n"))
  
  #Cannot check whether MLE has been reached for Poisson
  if(x$method == "poisson"){
    cat("Cannot check whether NPMLE has been reached for 'method = poisson'.")
  } else{
    mle_reached <- ifelse(x$isMLE, "", "NOT ")
    cat(paste0("NPMLE has ", mle_reached, "been reached. "))
    if(!x$isMLE){
      cat("Consider reducing 'tol' or increasing 'checkMLE_tol'.\n")
    } else{
      cat("\n")
    }
    cat(paste0("Conclusion based on checking reduced_gradient < ", x$checkMLE_tol))
  }
  
  
}