#' Sample from a markov chain multi state model with exactly observed transition times
#' 
#' @description Given a markov chain multi state model with exactly observed transition times,
#' sample from this chain at the observation times, giving interval censored observations (panel data).
#' 
#' 
#' @param time Times at which a transition occurs
#' @param stepf States at which the chain is in at \code{time}s
#' @param newtime Observation times of the chain, to create observed states
#' @param subst State to return if observation time is before first transition time. Default = -Inf.
#' @param to.data.frame Should the result be returned as a \code{data.frame}?
#' 
#' @author Hein Putter
#' 
#' @return A numeric \code{vector} or \code{data.frame} (if \code{to.data.frame = TRUE})
#' containing either the observed states or  the named columns \code{newtime} and 
#' \code{res}, representing the observation times and observed states.
#' 
#' @export
#' 
#' @examples
#' obs_states <- evalstep(time = seq(0, 20, 2), stepf = sample(1:9, 11, replace = TRUE),
#'                 newtime = c(-1, 1, 7, 9, 19))
#' obs_states


evalstep <- function (time, stepf, newtime, subst = -Inf, to.data.frame = FALSE) 
{
  #########################Parameter definitions################################
  n <- length(time)
  
  ########################Argument checks#######################################
  if (is.vector(stepf)) 
    if (length(stepf) != n) 
      stop("arguments 'time' and 'stepf' should have the same length")
  if (is.matrix(stepf) | is.data.frame(stepf)) 
    if (nrow(stepf) != n) 
      stop("argument 'stepf' should have the same number of rows as length of argument 'time'")
  if (any(!(order(time) == 1:n))) 
    stop("argument 'time' should be ordered")
  if (any(duplicated(time))) 
    stop("argument 'time' should not contain duplicates")
  if (any(is.infinite(time))) 
    stop("(-) infinity not allowed in 'time'")
  
  ########################Function body#########################################
  #Determine in which state the chain is at the observation times by cutting the 
  #transition time frame into pieces
  idx <- cut(newtime, c(-Inf, time, Inf), right = FALSE)
  idx <- as.numeric(idx)
  if (is.vector(stepf)) 
    res <- c(subst, stepf)[idx]
  if (is.matrix(stepf) | is.data.frame(stepf)) {
    stepf <- rbind(rep(subst, ncol(stepf)), stepf)
    res <- stepf[idx, ]
  }
  if (to.data.frame) 
    return(data.frame(newtime = newtime, res = res))
  else return(res)
}