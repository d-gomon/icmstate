#' Wrapper for the \code{\link[mstate:probtrans]{probtrans}} function
#' 
#' @description For \code{'smoothmsm'} objects: Determine transition probabilities for an \code{'smoothmsm'} object 
#' using the \code{\link[mstate:probtrans]{probtrans}} function.
#' 
#' @rdname transprob
#' 
#' @param object An "smoothmsm" object 
#' @param ... Further arguments to \code{\link[mstate:probtrans]{probtrans}}
#' 
#' @importFrom mstate probtrans
#' @importFrom utils modifyList
#' @method transprob smoothmsm
#' @export
#' 
#' 
#' 




transprob.smoothmsm <- function(object, ...){
  dots <- list(...)
  
  if ("variance" %in% names(dots)) {
    dots$variance <- FALSE
  } else {
    dots <- modifyList(dots, list(variance = FALSE))
  }
  do.call(mstate::probtrans, c(list(object = object$smoothmsfit), dots))
}