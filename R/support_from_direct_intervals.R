#' Estimate support of multiple transitions given direct transition intervals
#' 
#' @description Given only direct transition intervals, determine the support
#' for each transition separately using Hudgens(2001) result. Each state is
#' considered from a competing-risks viewpoint. Hudgens(2005) result is applied 
#' to see if the NPMLE for any of the transitions does not exist.
#' 
#' 
#' @param direct_intervals Output from \code{\link{direct_from_observed_intervals}}.
#' @param tmat A transition matrix as created by \code{transMat}
#'
#' @return A list containing the estimated support sets for each possible transition
#' in \code{tmat}.
#' 
#' 
#' @importFrom mstate to.trans2
#' 
#' @export
#' 
#' 




support_from_direct_intervals <- function(direct_intervals, tmat){
  
  # Remove CRAN notes
  from <- to <- NULL
  
  #For each possible transition, we want to apply Hudgens(2005) (implemented in 
  #supportHudgens()) on the determined direct transition intervals.
  
  #Additional properties extraction
  tmat2 <- to.trans2(tmat)
  M <- nrow(tmat2)
  
  #Initialize list to store supports
  transition_support_storage <- vector(mode = "list", length = M)
  
  
  #For each transition separately
  for(m in 1:M){
    state_from <- tmat2$from[m]
    state_to <- tmat2$to[m]
    
    #Extract only the necessary direct intervals for current transition.
    temp_intervals <- subset(direct_intervals, from == state_from & to == state_to)
    #If transition is not observed, skip to next.
    if(nrow(temp_intervals) == 0){
      transition_support_storage[[m]] <- list(graph = NULL, support = NULL, dir_graph = NULL, exist_MLE = FALSE)
      warning(paste0("No transitions of type ", m, " observed. Unclear if NPMLE exists."))
      next
    }
    
    #We can get multiple possible paths from a single observed transition
    which_duplicated <- which(duplicated(temp_intervals$id))
    if(length(which_duplicated) > 0){
      temp_intervals$id[which_duplicated] <- max(temp_intervals$id) + 1:length(which_duplicated)
    }
    
    #Initialize data.frame to use in supportHudgens()
    intervals <- matrix(NA, nrow = nrow(temp_intervals) * 2, ncol = 4)
    #browser()
    #For each interval, create a left-truncation intervals and a censoring interval
    for(j in 1:nrow(temp_intervals)){
      current_transition <- temp_intervals[j,]
      intervals[2*(j-1) + 1,] <- c(current_transition$time_from,
                                   current_transition$time_to, 
                                   1,
                                   current_transition$id)
      intervals[2*(j-1) + 2,] <- c(current_transition$entry_time,
                                   Inf,
                                   0,
                                   current_transition$id)
      
    }
    intervals <- as.data.frame(intervals)
    colnames(intervals) <- c("L", "R", "cens", "id")
    
    #Remove duplicated (can happen sometimes)
    #intervals <- intervals[!duplicated(intervals),]
    
    
    transition_support_storage[[m]] <- supportHudgens(intervals, reduction = TRUE, existence = TRUE)
    if(!transition_support_storage[[m]]$exist_mle){
      warning(paste0("NPMLE for transition ", m, " may not exist. Please check support set 
                     when running npmsm() without estimating support and compare the likelihood values."))
    }
  }
  
  
  return(transition_support_storage)
}























