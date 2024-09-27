#' @title Plot a "npmsm" object
#' 
#' @description A wrapper for \code{\link[mstate:msfit]{msfit}}
#' 
#' @param x An object of class \code{"npmsm"}
#' @param ... Additional arguments to \code{\link[mstate:msfit]{msfit}}
#' 
#' @return A plot will be produced in the plotting window
#'
#' @export
#' @import mstate
#'



plot.npmsm <- function(x, ...){
  plot(x$A, ...)
}