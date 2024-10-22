#' Get transition intervals from specified data
#' 
#' @description Given a sample from a multi-state model, summarize the transitions 
#' that have been observed.
#' 
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param tmat A transition matrix as created by \code{transMat}
#' 
#' @return A \code{data.frame} with the following named columns
#' \describe{
#'   \item{\code{entry_time}:}{Time of entry into "from" state;}
#'   \item{\code{time_from}:}{Last time subject(id) was seen in state "from";}
#'   \item{\code{time_to}:}{First time subject(id) was seen in state "to";}
#'   \item{\code{from}:}{State from which a transition was observed;}
#'   \item{\code{to}:}{State to which the transition was observed;}
#'   \item{\code{id}:}{Subject identifier;}
#' } For right-censored observations, entry_time denotes the first time
#' seen in the censored state, time_from the last time seen in the censored state,
#' time_to is \code{Inf}, from the censored state and to is \code{NA}.
#' 
#' 
#' 
#' @export
#' 
#' 




get_trans_intervals <- function(gd, tmat){
  
  # Remove CRAN check notes
  id <- NULL
  
  #Sort gd according to id and then time (neater this way)
  gd <- gd[order(gd["id"], gd["time"]),]
  
  
  #Extract some characteristics of the data
  unique_id <- unique(gd$id)
  n <- length(unique_id)
  #If tmat is supplied, determine which states are absorbing
  if(!missing(tmat)){
    tmat_temp <- tmat
    diag(tmat_temp) <- NA
    absorbing_states <- which(apply(tmat_temp, 1, function(x) all(is.na(x))))  
  }
  transient_states <- setdiff(1:nrow(tmat), absorbing_states)
  
  
  out <- NULL
  unique_ids <- unique(gd$id)
  for(i in unique_ids){
    gdi <- subset(gd, id == i)
    #Did we go into a new state in corresponding row? 
    new_state_times_idx <- which(c(FALSE, diff(gdi$state) != 0))
    new_state_times <- gdi$time[new_state_times_idx]
    
    #Get data to add to output df
    entry_time <- c(0, new_state_times[-length(new_state_times)])
    if(length(new_state_times) > 0){
      new_states <- gdi$state[new_state_times_idx]
      old_states <- gdi$state[new_state_times_idx - 1]
      #Add transitions to output data.frame
      out <- rbind(out, matrix(c(entry_time, gdi$time[new_state_times_idx - 1], gdi$time[new_state_times_idx], old_states, new_states, rep(i, length(old_states))), ncol = 6))
      #CHECK FOR RIGHT CENSORING
      final_state <- gdi$state[nrow(gdi)]
      last_entered_state <- new_states[length(new_states)]
      #If we are right-censored in a transient state
      if(final_state == last_entered_state & final_state %in% transient_states){
        final_state_first_time <- gdi$time[new_state_times_idx[length(new_state_times_idx)]]
        final_state_last_time <- gdi$time[nrow(gdi)]
        out <- rbind(out, c(final_state_first_time, final_state_last_time, Inf, final_state, NA, i))
      }
    } else{
      initial_state <- gdi$state[1]
      #If no transitions, must be censored in initial state
      if(initial_state %in% transient_states){
        out <- unname(rbind(out, c(gdi$time[1], gdi$time[nrow(gdi)], Inf, initial_state, NA, i)))
      } #Notice that there is no else, because if we have observed someone
      #initially in an absorbing state, this does not contribute to support sets (or even the likelihood)
    }
  }
  out <- as.data.frame(out)
  colnames(out) <- c("entry_time", "time_from", "time_to", "from", "to", "id")
  return(out)
}














