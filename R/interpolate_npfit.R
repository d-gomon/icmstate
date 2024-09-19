#' Given a \code{msfit} object, linearly interpolate the cumulative hazard 
#' taking into account the support sets for \code{msfit} objects.
#'
#' @description After using this function, use probtrans to get interpolated 
#' transition probabilities.
#'
#' @param msfit A \code{msfit} object.
#' @param times Times at which to interpolate the \code{msfit} object.
#'
#' @importFrom stats approx
#'
#' @return An \code{msfit} object containing the interpolated hazards
#' @export
#' 
#' 
#' @examples 
#' library(mstate)
#' tmat <- trans.illdeath()
#' times <- seq(0, 5, 0.1)
#' ms_fit <- list(Haz = data.frame(time = rep(times, 3),
#'                                 Haz = c(replicate(3, cumsum(runif(length(times), 0, 0.02)))),
#'                                 trans = rep(1:3, each = length(times))),
#'                trans = tmat)
#' class(ms_fit) <- "msfit"
#' 
#' ms_fit_interpolated <- interpol_msfit(ms_fit, seq(0, 5, 0.01))
#'


interpol_msfit <- function(msfit, times){
  numeric_support <- support_npmsm(msfit)
  #Number of possible transitions
  M <- length(numeric_support)
  out <- NULL
  #We interpolate for each transition separately.
  for (i in 1:M){
    Haz <- rbind(c(0, 0, i), msfit$Haz[msfit$Haz$trans == i,])
    supp <- numeric_support[[i]]$support[, c(1,2), drop = FALSE]
    
    #If there is no support set
    if(nrow(supp) == 0){
      supp <- rbind(supp, c(min(Haz$time), max(Haz$time)))
    }
    
    #Make interpolation grid
    x <- matrix(t(supp), nrow = 1, byrow = FALSE)
    y <-  Haz$Haz[findInterval(x, Haz$time)]
    x <- c(0, x)
    y <- c(0, y)
    
    out <- suppressWarnings({
      rbind(out, data.frame(time = times, Haz = stats::approx(x, y, xout = times, rule = 2)[["y"]], trans = rep(i, length(times)) ))
    })
      
  }
  
  out <- list(Haz = out,
              trans = msfit$trans)
  class(out) <- "msfit"
  return(out)
}















