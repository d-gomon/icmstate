#' @title Print a "npmsm" object
#' 
#' @description Print some details of a \code{\link{npmsm}} fit
#' 
#' @param object An object of class \code{"npmsm"}
#' @param ... Additional arguments to print
#' 
#' @return A summary of the fitted model will be displayed in the console
#'
#' @export
#'



print.npmsm <- function(object, ...){
  method_used <- ifelse(object$method == "multinomial", "Multinomial", "Poisson")
  crit_used <- switch(object$conv_crit, 
                      "haz" = "change in estimated intensities",
                      "prob" = "change in estimated probabilities",
                      "lik" = "change in observed-data likelihood")
  cat(paste0(method_used, " complete-data likelihood maximized through EM algorithm.\n"))
  cat(paste0("Convergence criterion: ", crit_used, " < ", object$tol, "\n"))
  if(object$it < object$maxit){
    cat(paste0("Algorithm converged after ", object$it, " out of max ", object$maxit, " iterations.\n"))  
  } else{
    cat(paste0("Algorithm failed to converge after ", min(object$it, object$maxit), " out of max ", object$maxit, " iterations.\n"))
    cat("Consider increasing 'maxit'.\n")
  }
  cat(paste0("Log likelihood value: ", object$ll, "\n"))
  
  #Cannot check whether MLE has been reached for Poisson
  if(object$method == "poisson"){
    cat("Cannot check whether NPMLE has been reached for 'method = poisson'.")
  } else{
    mle_reached <- ifelse(object$isMLE, "", "NOT ")
    cat(paste0("NPMLE has ", mle_reached, "been reached. "))
    if(!object$isMLE){
      cat("Consider reducing 'tol' or increasing 'checkMLE_tol'.\n")
    } else{
      cat("\n")
    }
    cat(paste0("Conclusion based on checking reduced_gradient < ", object$checkMLE_tol))
  }
  
  
}