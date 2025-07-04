#' Translate observed transition intervals into direct transition intervals
#' 
#' @description Given observed transition intervals, determine the "worst" (least informative) possible 
#' direct transition intervals that could have occurred to form this sample.
#' 
#' 
#' @param observed_intervals Output from \code{\link{get_trans_intervals}}.
#' @param tmat A transition matrix as created by \code{transMat}
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
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
#' @importFrom mstate to.trans2
#' @importFrom igraph graph_from_adjacency_matrix all_simple_paths
#' @export
#' 
#' 



direct_from_observed_intervals <- function(observed_intervals, tmat, gd){
  
  #Additional properties extraction
  tmat2 <- to.trans2(tmat)
  n_obs_trans <- nrow(observed_intervals)
  M <- nrow(tmat2)
  
  taus <- sort(unique(gd$time))
  K <- length(taus)
  
  #Initiate data frame to store intervals
  out_df <- NULL
  #out_df must be in a format usable by supportHudgens()
  #columns: L, R, cens, id
  
  #Create graph to check possible transitions
  adjacency_matrix <- tmat
  adjacency_matrix[!is.na(tmat)] <- 1
  adjacency_matrix[is.na(tmat)] <- 0
  colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- 1:nrow(adjacency_matrix)
  adjacency_graph <- graph_from_adjacency_matrix(adjacency_matrix)
  
  #Iterate over each observed transition
  for(i in 1:n_obs_trans){
    #We look at a single transition at a time
    current_trans <- observed_intervals[i,]
    #columns: entry_time, time_from, time_to, from, to, id
    
    #We first cover the right-censored case
    if(is.na(current_trans$to)){
      possible_future_states <- tmat2$to[which(tmat2$from == current_trans$from)]
      for(st in possible_future_states){
        out_df <- rbind(out_df, c(current_trans$entry_time,
                                  current_trans$time_from,
                                  current_trans$time_to,
                                  current_trans$from,
                                  st,
                                  current_trans$id))
      }
      next
    }
    
    
    #If we do not want direct transitions to contribute to indirect
    #
    #
    # #We check if current_transition could have been taken directly
    # direct_transition <- !is.na(adjacency_matrix[current_trans$from, current_trans$to])
    # 
    # #If direct transition is possible, consider only that transition
    # if(direct_transition){
    #   possible_paths <- list(c(current_trans$from, current_trans$to))
    # } else{
    #   #Determine all possible paths this transition could have occurred in:
    #   possible_paths <- all_simple_paths(adjacency_graph, from = current_trans$from,
    #                                      to = current_trans$to)
    # }
    
    possible_paths <- all_simple_paths(adjacency_graph, from = current_trans$from,
                                       to = current_trans$to)
    
    #Get indices corresponding to taus at which the transition has occurred
    tau_idx_time_from <- which(taus == current_trans$time_from)
    tau_idx_time_to <- which(taus == current_trans$time_to)
    
    for(j in 1:length(possible_paths)){
      #Consider one of the possible paths
      path <- possible_paths[[j]]
      len_path <- length(path)
      #Iterate over all direct transition within path
      for(k in 1:(len_path-1)){
        #Consider a single direct transition in this path
        trans <- path[c(k, k+1)]
        #Number of transitions to the left and right of current transition
        n_trans_left <- k-1
        n_trans_right <- (len_path-1)-k
        
        #When there are not enough bins to incorporate this, we need to adjust
        #to least possible information
        if(n_trans_left + n_trans_right >= (tau_idx_time_to - tau_idx_time_from) -1 ){
          n_trans_left <- n_trans_right <- 0
        }
        
        #Find new time indices for "worst case" transitions
        tau_idx_left <- tau_idx_time_from + n_trans_left
        tau_idx_right <- tau_idx_time_to - n_trans_right
        
        #First transition is left-truncated at initial entry time.
        if(k == 1){ 
          #Add intervals to output
          #entry_time time_from time_to from to id
          out_df <- rbind(out_df, c(current_trans$entry_time,
                                    taus[tau_idx_left],
                                    taus[tau_idx_right],
                                    trans[1],
                                    trans[2],
                                    current_trans$id))
        } else{ #Other transitions are "worst-case" left-truncated
          out_df <- rbind(out_df, c(taus[tau_idx_left],
                                    taus[tau_idx_left],
                                    taus[tau_idx_right],
                                    trans[1],
                                    trans[2],
                                    current_trans$id))
        }
      }
    }
  }
  out_df <- as.data.frame(out_df)
  colnames(out_df) <- c("entry_time", "time_from", "time_to", "from", "to", "id")
  
  return(out_df)
}










