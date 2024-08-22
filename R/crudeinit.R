#' Find crude initial estimates for the transition intensities in an interval-
#' censored Markov Multi State model without loops.
#' 
#' 
#' @description The crude initial estimates are determined by dividing the 
#' number of people 
#' 
#' 
#' @param gd A \code{data.frame} with the following named columns
#'\describe{
#'   \item{\code{id}:}{Subject idenitifier;}
#'   \item{\code{state}:}{State at which the subject is observed at \code{time};}
#'   \item{\code{time}:}{Time at which the subject is observed;}
#' } The true transition time between states is then interval censored between the times.
#' @param gd_original The original \code{gd} data.frame, only used for exactly 
#' observed event times.
#' @param tmat A transition matrix as created by \code{transMat}
#' 
#' @importFrom mstate to.trans2
#' @importFrom igraph graph_from_adjacency_matrix all_simple_paths
#' @importFrom prodlim row.match
#' 
#' 
#' @keywords internal
#' @noRd
#' 


crudeinit.npmsm <- function(gd, taus, tmat){
  #gd must be gd_original in the exactly observed case!
  if(missing(taus)){
    taus <- sort(unique(gd[, "time"]))
    taus <- taus[-1] #Minimal time of gd has been set to 0
  }
  
  tmat2 <- to.trans2(tmat)
  
  direct_transitions <- tmat2[, c(2, 3)]
  
  
  #Create graph to check possible transitions
  adjacency_matrix <- tmat
  adjacency_matrix[!is.na(tmat)] <- 1
  colnames(adjacency_matrix) <- rownames(adjacency_matrix) <- 1:nrow(adjacency_matrix)
  adjacency_graph <- graph_from_adjacency_matrix(adjacency_matrix)
  
  #Get all transition intervals
  trans_intervals <- get_trans_intervals(gd, tmat)
  #Determine the unique transition intervals
  unique_transitions <- unique(trans_intervals[, c("from", "to")])
  
  #Determine to which direct transitions the unique transition intervals contribute
  trans_contribution <- vector(mode = "list", length = nrow(unique_transitions))
  for(i in 1:nrow(unique_transitions)){
    current_unique_transition <- unique_transitions[i,]
    if(is.na(current_unique_transition[, "to"])){ #If right-censored, no contribution to transitions
      trans_contribution[[i]] <- -1
    } else{
      possible_paths <- all_simple_paths(adjacency_graph, from = current_unique_transition[, "from"],
                                         to = current_unique_transition[, "to"])
      #Check for each path which direct transitions are contained therein
      for(j in 1:length(possible_paths)){
        current_path <- possible_paths[[j]]
        for(k in 1:(length(current_path)-1))
        trans_contribution[[i]] <- c(trans_contribution[[i]], prodlim::row.match(current_path[c(k, k+1)], direct_transitions))
      }
    }
    trans_contribution[[i]] <- sort(trans_contribution[[i]])
  }
  
  
  K <- length(taus)
  M <- nrow(tmat2)
  n <- length(unique(gd[, "id"]))
  
  d_init <- array(0, dim = c(K, M, n))
  Y_init <- array(0, dim = c(K, M, n))
  
  
  for(i in 1:nrow(trans_intervals)){
    #Current transition for a single subject
    current_trans <- trans_intervals[i,]
    #Which times does it contribute to?
    wh <- which(taus > current_trans[, "time_from"] & taus <= current_trans[, "time_to"])
    
    #Which direct transitions are possible?
    corresponding_unique_transition <- prodlim::row.match(current_trans[, c("from", "to")], unique_transitions)
    contribute_transitions <- trans_contribution[[corresponding_unique_transition]]
    
    d_init[wh, contribute_transition, current_trans[, "id"]] <- 1
    
    
  }
  
}







