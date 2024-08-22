#' Given a \code{msm} object, determine the transition probabilities (forward) 
#' at the specified times, starting from the smallest specified time.
#'
#' @description The output is similar to probtrans from the mstate package.
#'
#' @param msm A \code{msm} object.
#' @param trans A transition number for which the interpolated hazard should be 
#' returned
#'
#' @return A \code{probtrans} object containing the estimated transition probabilities.
#'
#' @export
#' @noRd
#'


probtrans_msm <- function(msm, times){
  #Number of allowed transitions
  M <- msm$qmodel$npars
  #Number of states
  n_states <- msm$qmodel$nstates
  
  #We want to get output as from probtrans
  #For this we apply the pmatrix.msm function many times
  times <- times - min(times)
  out <- lapply(times, function(y) msm::pmatrix.msm(x = msm, t = y, t1 = 0))
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
  res$trans <- tmat_from_msm(msm)
  return(res)
}