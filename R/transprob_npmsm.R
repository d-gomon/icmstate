#' Wrapper for the \code{\link[mstate:probtrans]{probtrans}} function
#' 
#' @description For \code{'npmsm'} objects: Determine transition probabilities for an \code{'npmsm'} object 
#' using the \code{\link[mstate:probtrans]{probtrans}} function.
#' 
#' @rdname transprob
#' 
#' @param object An "npmsm" object 
#' @param ... Further arguments to \code{\link[mstate:probtrans]{probtrans}}
#' 
#' @importFrom mstate probtrans
#' @importFrom utils modifyList
#' @method transprob npmsm
#' @export
#' 
#' 
#' 




transprob.npmsm <- function(object, ...){
  dots <- list(...)
  
  if ("variance" %in% names(dots)) {
    dots$variance <- FALSE
  } else {
    dots <- modifyList(dots, list(variance = FALSE))
  }
  do.call(mstate::probtrans, c(list(object = object$A), dots))
}



#' Wrapper for \code{\link[mstate:probtrans]{probtrans}}
#' 
#' 
#' @rdname transprob
#' @description
#' For \code{'msfit'} objects: Wrapper for \code{\link[mstate:probtrans]{probtrans}}
#' 
#' 
#' @importFrom mstate probtrans
#' 
#' @param object Object of class 'msfit'
#' @param ... Further arguments to \code{\link[mstate:probtrans]{probtrans}}
#' 
#' @method transprob msfit
#' @export
#' 
transprob.msfit <- function(object, ...){
  mstate::probtrans(object, ...)
}


#' S3 method for 'transprob'
#' 
#' @details
#' Can be used for objects of class 'npmsm', 'msm' and 'msfit'
#' 
#' @param object Object of compatible class
#' @param ... Further arguments to \code{\link[mstate:probtrans]{probtrans}}
#' 
#' 
#' @export
#'
#'

transprob <- function(object, ...){
  UseMethod("transprob")
}