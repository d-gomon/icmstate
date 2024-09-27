#'
#' @description For \code{'msm'} objects: determine transition probabilities
#' (as in \code{\link[mstate:probtrans]{probtrans}})  from an 
#' \code{\link[msm:msm]{msm}} object.
#' 
#' @rdname transprob
#'
#' @param object A \code{msm} object.
#' @param predt A positive number indicating the prediction time. This is 
#'  the time at which the prediction is made. If missing, smallest time of 
#' \code{times} is chosen.
#' @param times A vector of times at which the transition probabilities should 
#' be determined.
#' @param ... Further arguments to transprob
#'
#' @return A \code{probtrans} object containing the estimated transition probabilities.
#' @method transprob msm
#' @export
#' @importFrom msm pmatrix.msm
#'


transprob.msm <- function(object, predt, times, ...){
  if(missing(predt)){
    predt = min(times)
  }
  #Number of allowed transitions
  M <- object$qmodel$npars
  #Number of states
  n_states <- object$qmodel$nstates
  
  #We want to get output as from probtrans
  #For this we apply the pmatrix.msm function many times
  times <- times - min(times)
  out <- lapply(times, function(y) msm::pmatrix.msm(x = object, t = y, t1 = predt))
  #Faster if you just give times to t argument, but the output is some weird object that I don't know how to manipulate
  
  #For now we store our data in an array, to later transform it into a list as in probtrans
  out_arr <- array(NA, dim = c(n_states, n_states + 1, length(times)))
  
  for(i in 1:length(out)){
    out_arr[, 1, i] <- times[i]
    out_arr[, 2:(n_states + 1) , i] <- out[[i]]
  }
  res <- vector(mode = "list", length = n_states)
  for(s in 1:n_states){
    res[[s]] <- as.data.frame(t(out_arr[s, , ]))
    colnames(res[[s]]) <- c("time", paste("pstate", 1:n_states, sep = ""))
  } 
  res$trans <- tmat_from_msm(object)
  class(res) <- "probtrans"
  return(res)
}