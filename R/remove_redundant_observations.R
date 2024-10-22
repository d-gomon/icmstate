#' Remove redundant observations from supplied data frame
#' 
#' @description Remove redundant observed states from a supplied data frame. 
#' Observations are redundant either when we observe an absorbing state multiple 
#' times (as we cannot leave an absorbing state), or 
#' when a transient state is observed multiple times between transitions (as 
#' we cannot have loops, therefore no extra information is provided when 
#' we observe a transient state multiple times).
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
#' @export
#' 
#' 
#' @returns A \code{data.frame} containing the information contained in  the 
#' input \code{data.frame} \code{gd},
#' but without redundant observations. Depending on whether \code{tmat} was 
#' specified the function may remove more observations.
#' 
#' @importFrom utils tail
#' 
#' @examples
#' #We simulate some data
#' #Function to generate evaluation times: at 0 and uniform inter-observation
#' eval_times <- function(n_obs, stop_time){
#'  cumsum( c( 0,  runif( n_obs-1, 0, 2*(stop_time-4)/(n_obs-1) ) ) )
#' }
#' 
#' #Simulate illness-death model data with Weibull transitions
#' sim_dat <- sim_id_weib(n = 20, n_obs = 6, stop_time = 15, eval_times = eval_times,
#'                        start_state = "stable", shape = c(0.5, 0.5, 2), 
#'                        scale = c(5, 10, 10/gamma(1.5)))
#' visualise_msm(sim_dat)
#' require(mstate)
#' sim_dat_clean <- remove_redundant_observations(sim_dat, trans.illdeath())
#' visualise_msm(sim_dat_clean)
#' 




remove_redundant_observations <- function(gd, tmat){
  
  #Sort gd according to id and then time (neater this way)
  gd <- gd[order(gd[, "id"], gd[, "time"]),]
  
  #Extract some characteristics of the data
  unique_id <- unique(gd[, "id"])
  n <- length(unique_id)
  #If tmat is supplied, determine which states are absorbing
  if(!missing(tmat)){
    tmat_temp <- tmat
    diag(tmat_temp) <- NA
    absorbing_states <- which(apply(tmat_temp, 1, function(x) all(is.na(x))))  
  }
  
  which_duplicated <- which(duplicated(gd))
  if(anyDuplicated(gd)){
    which_duplicated <- which(duplicated(gd))
    warning(paste0("Duplicated entries found for id's: ", toString(unique(gd[which_duplicated, "id"])), ". 
                   Removing duplicated, but consider checking the data processing procedure."))
    gd <- gd[-which_duplicated, ]
  }
  
  # Remove duplicated observations from gd when:
  # - Landed in an absorbing state (only if tmat is present)
  # - Observed intermediate state multiple times (we cannot have loops, so must have stayed in this state)
  

  #Keep track of indices to remove
  idx_remove <- NULL
  #Loop over subjects  
  for(i in 1:n){
    i_idx <- which(gd[, "id"] == unique_id[i])
    smallest_entry <- min(i_idx)
    #Which state entries are duplicated for current subject?
    bool_dupl <- duplicated(gd[i_idx, "state"])
    duplicated_idx <- bool_dupl
    #We want to retain the last duplicated entry of each state,
    #because we transition to a non-duplicated state afterwards
    duplicated_idx[which(diff(bool_dupl) < 0)] <- FALSE
    #We also retain the final entry if it is a non-absorbing state, as this gives us information
    if(!missing(tmat)){
      if(!(gd[tail(i_idx, 1), "state"] %in% absorbing_states)){
        duplicated_idx[length(duplicated_idx)] <- FALSE
      }      
    } else{ #If we do not know whether we have absorbing states, we also retain the final state
      duplicated_idx[length(duplicated_idx)] <- FALSE
    }
    #Remove entries 
    remove_idx <- which(duplicated_idx)
    if(length(remove_idx) > 0){
      idx_remove <- c(idx_remove, smallest_entry-1 + remove_idx)
      #gd <- gd[-i_idx[remove_idx],]
    }
  }
  if(length(idx_remove) > 0){
    gd <- gd[-idx_remove,]  
  }
  return(gd)
}