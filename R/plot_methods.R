#' @title Plot a "npmsm" object
#' 
#' @description Plot the cumulative intensities of a \code{'npmsm'} objects.
#' A wrapper for \code{\link[mstate:plot.msfit]{plot.msfit}} from the
#' \code{\link[mstate:mstate]{mstate}} package. 
#' 
#' @param x An object of class \code{"npmsm"}
#' @param ... Additional arguments to \code{\link[mstate:msfit]{msfit}}
#' 
#' @return A plot will be produced in the plotting window
#'
#' @method plot npmsm
#'
#' @export
#' @import mstate
#'



plot.npmsm <- function(x, ...){
  plot(x$A, ...)
}


#' @title Plot a "smoothmsm" object
#' 
#' @description Plot the cumulative intensities of a \code{'smoothmsm'} objects.
#' For the cumulative hazard, a wrapper for \code{\link[mstate:plot.msfit]{plot.msfit}} from the
#' \code{\link[mstate:mstate]{mstate}} package. 
#' 
#' @param x An object of class \code{"smoothmsm"}
#' @param haz Type of hazard to plot. Choices are \code{"haz"} (default: hazard),
#' \code{"cumhaz"} (cumulative hazard), \code{"loghaz"} (log-hazard).
#' @param ... Additional arguments to \code{\link[mstate:plot.msfit]{plot.msfit}}
#' 
#' @return A plot will be produced in the plotting window
#'
#' @method plot smoothmsm
#' @import ggplot2
#'
#' @export
#'



plot.smoothmsm <- function(x, haz = c("haz", "cumhaz", "loghaz"), ...){
  haz <- match.arg(haz, choices = c("haz", "cumhaz", "loghaz"))
  #Everything is basically already pre-calculated in create_smoothmsfit
  msfit_obj <- x$smoothmsfit
  args <- list(...)
  if(haz == "haz"){
    #Plot haz
    msfit_obj$Haz <- msfit_obj$haz
    if(!"ylab" %in% names(args)){
      ylab <- "hazard"
      args <- c(args, ylab = ylab)
    }
  } else if(haz == "cumhaz"){
    #Do nothing, because Haz is the standard
    if(!"ylab" %in% names(args)){
      ylab <- "Cumulative hazard"
      args <- c(args, ylab = ylab)
    }
  } else{
    #Plot log-haz
    msfit_obj$Haz <- msfit_obj$loghaz
    if(!"ylim" %in% names(args)){
      ylim = c(min(msfit_obj$Haz$haz), max(msfit_obj$Haz$haz))
      args <- c(args, ylim = list(ylim))
    }
    if(!"ylab" %in% names(args)){
      ylab <- "log-hazard"
      args <- c(args, ylab = ylab)
    }
  }
  args <- c(x = list(msfit_obj), args)
  do.call(plot, args)
}







