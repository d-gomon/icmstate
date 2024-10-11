#' Simulate panel data from an illness-death model with weibull transitions
#' 
#' @description An illness-death model has 3 transitions:
#' \describe{
#' \item{\code{1}:}{State 1 (Healthy) to State 2 (Illness);}
#' \item{\code{2}:}{State 1 (Healthy) to State 3 (Death);}
#' \item{\code{3}:}{State 2 (Illness) to State 3 (Death);}
#' }
#' Using this function, it is possible to simulate data from an illness-death 
#' model with Weibull transition intensities. Requires the use of an external 
#' (self-written) function to generate observation times.
#' 
#' @details Taking \code{shape = 1} we get an exponential distribution with rate
#' \code{1/scale}
#' 
#' @param n Number of subjects to generate paths for.
#' @param n_obs Number of observations in time period for each subject.
#' @param stop_time Largest time at which the model is considered.
#' @param eval_times A function which returns the evaluation times for a subject.
#' Must have as arguments at least \code{n_obs} and \code{stop_time}.
#' @param start_state In which states can subjects start? Either everyone starts 
#' in state 1 ("stable") or equal probability to start in state 1 or 2 ("equalprob").
#' @param shape Vector of shape parameters for the 3 transitions. See \code{\link[stats:rweibull]{rweibull}}.
#' The first entry will be used for the first transition and so on.
#' @param scale Vector of scale parameters for the 3 transitions. See \code{\link[stats:rweibull]{rweibull}}
#' The first entry will be used for the first transition and so on.
#' @param ... Further parameters to \code{eval_times} function.
#' 
#' @export
#' 
#' @importFrom stats rweibull runif
#' 
#' @examples 
#' #Function to generate evaluation times: at 0 and uniform inter-observation
#' eval_times <- function(n_obs, stop_time){
#'   cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
#' }
#' 
#' #Simulate illness-death model data with Weibull transitions
#' sim_dat <- sim_id_weib(n = 20, n_obs = 6, stop_time = 15, eval_times = eval_times,
#' start_state = "stable", shape = c(0.5, 0.5, 2), scale = c(5, 10, 10/gamma(1.5)))
#' 
#' visualise_msm(sim_dat)
#' 





sim_id_weib <- function(n, n_obs, stop_time, eval_times, start_state = c("stable", "equalprob"),
                        shape, scale, ...){
  out <- NULL
  start_state <- match.arg(start_state)
  for(i in 1:n){
    subj_eval_times <- eval_times(n_obs = n_obs, stop_time = stop_time, ...)
    subj_start_state <- ifelse(start_state == "stable", 1, sample(c(1, 2), size = 1, prob = c(0.5, 0.5)))
    subj_trans_times <- 0
    subj_states <- subj_start_state
    #If we start in 1, need to generate survival times to both state 2 and 3.
    if(subj_start_state == 1){
      times <- c(rweibull(1, shape = shape[1], scale = scale[1]), 
                 rweibull(1, shape = shape[2], scale = scale[2]))
      if(times[1] <= times[2]){
        #If we went to state 2 first, also need to generate survival time to state 3
        #However, we have entered state 2 at time times[1], so need to generate from that time on
        #For this we use the truncated inverse Weibull hazard defined below
        time_3 <- inverse_weib_haz(t = -log(runif(1, 0, 1)), s = times[1], shape = shape[3], scale = scale[3])
        
        subj_trans_times <- c(subj_trans_times, times[1], time_3)
        subj_states <- c(subj_states, 2, 3)
      } else if(times[2] < times[1]){ #Else only need to generate survival time to state 3.
        subj_trans_times <- c(subj_trans_times, times[2])
        subj_states <- c(subj_states, 3)
      }
    } else{
      time_3 <- rweibull(1, shape = shape[3], scale = scale[3])
      subj_trans_times <- c(subj_trans_times, time_3)
      subj_states <- c(subj_states, 3)
    } 
    observed_subj_states <- evalstep(time = subj_trans_times, stepf = subj_states,
                         newtime = subj_eval_times)
    out_subj <- matrix(NA, nrow = length(observed_subj_states), ncol = 3)
    out_subj[, 1] <- rep(i, length(observed_subj_states))
    out_subj[, 2] <- subj_eval_times
    out_subj[, 3] <- observed_subj_states
    out <- rbind(out, out_subj)
  }
  colnames(out) <- c("id", "time", "state")
  return(as.data.frame(out))
}



#' 
#' Inverse cumulative hazard at time t of left-truncated Weibull distribution at 
#' time s.
#' 
#' @description Usefull for simulating. Generate a sample from a truncated 
#' Weibull distribution at time s by taking t = -log(U) with U unif(0,1)
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' 

inverse_weib_haz <- function(t, s, shape, scale){
  scale * (t + (s/scale)^shape)^(1/shape)
}


#' 
#' Inverse cdf at time t of left-truncated Weibull distribution at 
#' time s.
#' 
#' @description Usefull for simulating. Generate a sample from a truncated 
#' Weibull distribution at time s by taking t = U with U unif(0,1)
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' 

inverse_weib_cdf <- function(t, s, shape, scale){
  scale * ((s/scale)^shape - log(1-t))^(1/shape)
}



#' sim_id_weib but sanity check to see if we can reconstruct the cumulative hazards
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' @importFrom stats rweibull runif
#' 
#' 
#' 


sim_id_weib_exact <- function(n, n_obs, stop_time, eval_times, start_state = c("stable", "equalprob"),
                        shape, scale, ...){
  out <- NULL
  start_state <- match.arg(start_state)
  for(i in 1:n){
    subj_eval_times <- eval_times(n_obs = n_obs, stop_time = stop_time, ...)
    subj_start_state <- ifelse(start_state == "stable", 1, sample(c(1, 2), size = 1, prob = c(0.5, 0.5)))
    subj_trans_times <- 0
    subj_states <- subj_start_state
    #If we start in 1, need to generate survival times to both state 2 and 3.
    if(subj_start_state == 1){
      times <- c(rweibull(1, shape = shape[1], scale = scale[1]), 
                 rweibull(1, shape = shape[2], scale = scale[2]))
      if(times[1] <= times[2]){
        #If we went to state 2 first, also need to generate survival time to state 3
        #However, we have entered state 2 at time times[1], so need to generate from that time on
        #For this we use the truncated inverse Weibull hazard defined below
        time_3 <- inverse_weib_haz(t = -log(runif(1, 0, 1)), s = times[1], shape = shape[3], scale = scale[3])
        #Or generate through inverse cdf, also a possibility
        #time_3 <- inverse_weib_cdf(t = runif(1, 0, 1), s = times[1], shape = shape[3], scale = scale[3])
        
        subj_trans_times <- c(subj_trans_times, times[1], time_3)
        subj_states <- c(subj_states, 2, 3)
      } else if(times[2] < times[1]){ #Else only need to generate survival time to state 3.
        subj_trans_times <- c(subj_trans_times, times[2])
        subj_states <- c(subj_states, 3)
      }
    } else{
      time_3 <- rweibull(1, shape = shape[3], scale = scale[3])
      subj_trans_times <- c(subj_trans_times, time_3)
      subj_states <- c(subj_states, 3)
    } 
    observed_subj_states <- evalstep(time = subj_trans_times, stepf = subj_states,
                                     newtime = subj_eval_times)
    out_subj <- matrix(NA, nrow = length(subj_trans_times), ncol = 3)
    out_subj[, 1] <- rep(i, length(subj_trans_times))
    out_subj[, 2] <- subj_trans_times
    out_subj[, 3] <- subj_states
    out <- rbind(out, out_subj)
  }
  colnames(out) <- c("id", "time", "state")
  return(as.data.frame(out))
}

