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
#' @description Plot the cumulative intensities of a \code{'npmsm'} objects.
#' A wrapper for \code{\link[mstate:plot.msfit]{plot.msfit}} from the
#' \code{\link[mstate:mstate]{mstate}} package. 
#' 
#' @param x An object of class \code{"smoothmsm"}
#' @param ... Additional arguments to \code{\link[base:plot]{plot}}
#' 
#' @return A plot will be produced in the plotting window
#'
#' @method plot smoothmsm
#' @import ggplot2
#'
#' @export
#'



plot.smoothmsm <- function(x, ...){
  plot(x$A, ...)
}






