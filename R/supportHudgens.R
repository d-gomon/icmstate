#' Determine the support of the NPMLE for interval censored data.
#' 
#' @description Given censoring/truncation intervals, find the maxcliques and determine
#'  the support of the interval censored problem.
#' 
#' @param intervals A data.frame with 3 columns containing half-open intervals (left open, right closed)
#' and an indicator whether the interval results from a censored transition
#'  or truncation:
#' \describe{
#'   \item{\code{L}:}{Left side of interval;}
#'   \item{\code{R}:}{Right side of interval;}
#'   \item{\code{cens}:}{Indicator whether interval resulted from interval censoring or left truncation
#'   (1 = censoring, 0 = truncation);}
#'   \item{\code{id}:}{(optional) Identifier for the observation this interval belongs to (numeric/integer).
#'    Only required if existence = TRUE;}
#' } Note that the truncation intervals need to be in the form (N, Inf] with N a numeric value.
#' @param reduction Should the support be reduced using Lemma 3 from Hudgens (2005)? This 
#' requires checking an extra condition. Default is TRUE. 
#' @param existence Should the existence of the NPMLE be checked using Theorem 1/Lemma 4 from 
#' Hudgens (2005)? Requires \code{id} to be present in \code{intervals}. Default is FALSE. 
#' 
#' @return 
#' \itemize{
#' \item \code{graph}: An \code{igraph} object representing the censoring/truncation intervals
#' \item \code{support}: Support estimated from the censoring intervals
#' \item \code{dir_graph}: A directed \code{igraph} object used to determine whether the NPMLE
#' exists in the presence of left-truncation.
#' \item \code{exist_mle}: Logical output indicating whether the NPMLE exists.
#' }
#' 
#' 
#' @references Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#' 
#' @import igraph
#' @export
#' 
#' 




supportHudgens <- function(intervals, reduction = TRUE, existence = FALSE){
  
  #Data checks are performed in helper functions
  if(!any(intervals$cens == 0)){
    #If no intervals are truncated, we can skip the reduction step!
    reduction = FALSE
    if(isTRUE(existence)){
      existence = FALSE
      message("Existence cannot be checked if there are no truncation intervals. (The NPMLE will exist, but might not be unique)")
    }
  }
  
  #We can only check existence if the reduction step has already been employed
  if(isTRUE(existence)){
    reduction = TRUE
  }
  
  #Construct a graph from the intervals according to Hudgens (2005)
  int_graph <- graphfromIntervals(intervals)
  
  #Find the maxcliques in the graph
  maxcliq <- max_cliques(int_graph)
  #Maxcliques is a list containing maxcliques
  #Each list element has elements $L, $R, $cens, $name
  
  
  
  #Determine which maxcliques contain a censoring interval:
  supp_intervals_idx <- which(lapply(maxcliq, function(x) sum(x$cens)) > 0)
  #Extract only the maxcliques containing censoring intervals
  maxcliq_censoring <- maxcliq[supp_intervals_idx]
  #How many support intervals do we have?
  n_supp_intervals <- length(supp_intervals_idx)
  #Initiate matrix to store support intervals
  supportdf <- matrix(NA, nrow = n_supp_intervals, ncol = 2)
  
  #Iterate over each maxclique to:
  #1. Determine if it contains a censoring interval
  #2. If 1 is true, then take the intersection of the intervals
  for(i in 1:n_supp_intervals){
    #Keep track of which maxclique is associated with support interval:
    j <- supp_intervals_idx[i]
    L_sup <- max(maxcliq_censoring[[i]]$L)
    R_sup <- min(maxcliq_censoring[[i]]$R)
    supportdf[i,] <- c(L_sup, R_sup)
  }
  
  #Lemma 3 Hudgens (2005)
  #Check whether we should and can reduce the support
  if(isTRUE(reduction)){

    #Find index of rightmost maxclique and store the maxclique
    rightmost_clique_idx <- which.max(supportdf[, 1])
    rightmost_clique <- maxcliq_censoring[[rightmost_clique_idx]]
    #Extract the censored (NOT TRUNCATED) vertices from this maxclique
    rightmost_clique_vertices <- rightmost_clique[which(rightmost_clique$cens == 1)]
    
    #Storage for removing indices:
    k_store <- NULL
    #Iterate over all maxcliques except the rightmost one.
    for(k in setdiff(1:n_supp_intervals, rightmost_clique_idx)){
      temp_clique <- maxcliq_censoring[[k]]
      temp_clique_vertices <- temp_clique[which(temp_clique$cens == 1)]
      if(all(temp_clique_vertices %in% rightmost_clique_vertices)){
        k_store <- c(k_store, k)
      }
    }
    #If we find a maxclique that meets Lemma 3 from Hudgens (2005)
    #we remove the corresponding interval from the support.
    if(!is.null(k_store)){
      supportdf <- supportdf[-k_store, , drop = FALSE]  
    }
  }
  
  #Post-processing of support output
  colnames(supportdf) = c("L", "R")
  supportdf <- supportdf[order(supportdf[,1]), , drop = FALSE]
  
  out <- list(graph = int_graph,
              support = supportdf)
  
  #Check Existence NPMLE
  #Hudgens (2005) Theorem 1 + Lemma 4
  if(isTRUE(existence)){
    exist <- existenceNPMLE(intervals, supportdf)
    out$dir_graph = exist$dir_graph
    out$exist_mle = exist$exist_mle
  }
  
    
  return(out)
  
}

#' Construct Graph from censoring/truncation intervals 
#' 
#' @description Given intervals, construct a graph containing vertices 
#' representing these intervals and edges between the vertices if the intervals 
#' intersect. See Hudgens (2005). 
#' 
#' @param intervals A data.frame with 3 columns containing half-open intervals (left open, right closed)
#' and an indicator whether the interval results from a censored transition
#'  or truncation:
#' #'\describe{
#'   \item{\code{L}:}{Left side of interval;}
#'   \item{\code{R}:}{Right side of interval;}
#'   \item{\code{cens}:}{Indicator whether interval resulted from censoring or truncation
#'   (1 = censoring, 0 = truncation);}
#' } Note that the truncation intervals need to be in the form (N, Inf] with N a numeric value.
#' 
#' 
#' @references Michael G. Hudgens, On Nonparametric Maximum Likelihood Estimation with 
#' Interval Censoring and Left Truncation, Journal of the Royal Statistical Society 
#' Series B: Statistical Methodology, Volume 67, Issue 4, September 2005, Pages 573-587,
#'  \doi{10.1111/j.1467-9868.2005.00516.x}
#' 
#' @import igraph
#' @export
graphfromIntervals <- function(intervals){
  #----------------Data checks-----------------------------------
  if(!all(intervals$R > intervals$L)){
    stop("Right endpoint of intervals must be larger than left endpoint.")
  }
  if(!all(c("L", "R", "cens") %in% colnames(intervals))){
    stop("Please provide a data.frame with at least columns named L, R and cens.")
  }
  intervals <- intervals[,c("L", "R", "cens")]
  
  
  
  #------------------Function body-------------------------------
  
  #If we have truncation: take complement of truncation interval in union of all truncation intervals
  if(any(intervals$cens == 0)){
    #Process truncation intervals into complements
    trunc_idx <- which(intervals$cens == 0)
    #Minimal value of left endpoint for truncation intervals
    Lmin_trunc <- min(intervals[trunc_idx,"L"])
    intervals[trunc_idx, "R"] <- intervals[trunc_idx, "L"]
    intervals[trunc_idx, "L"] <- Lmin_trunc
  }
  
  #Give the intervals a temporary name so we can match them later
  intervals["name"] <- 1:nrow(intervals)
  intervals <- intervals[, c("name", "L", "R", "cens")]
  
  #Two intervals (L1, R1] and (L2, R2] intersect if:
  #R2 > L1 and L2 < R1
  RgreaterL <- outer(intervals$R, intervals$L, ">")
  LsmallerR <- outer(intervals$L, intervals$R, "<")
  
  #Determine intersections (both have to be true, so just sum, if sum==2 then intersect)
  intersec <- (RgreaterL + LsmallerR)
  #Remove lower diagonal, as upper and lower are equivalent
  intersec[lower.tri(RgreaterL, diag = TRUE)] <- 1
  #We get the row and column indices to construct the graph
  rowcol_idx <- which(intersec == 2, arr.ind = TRUE)
  
  g <- graph_from_data_frame(rowcol_idx, directed = FALSE, vertices = intervals)
  #Name the vertices using Truncated/Censored: T/C and their L/R values
  V(g)$name <- apply(intervals, 1, function(x) paste0(ifelse(x[4] == 1, "C", "T"), round(x[2], 2), "-", round(x[3],2)))
  return(g)
}






#' Check existence of NPMLE
#' 
#' @description Using Theorem 1 from Hudgens (2005) we can check whether an 
#' NPMLE exists. This procedure is implemented in this function.
#' 
#' @param intervals A data.frame with 4 columns containing half-open intervals (left open, right closed)
#' and an indicator whether the interval results from a censored transition
#'  or truncation:
#' \describe{
#'   \item{\code{L}:}{Left side of interval;}
#'   \item{\code{R}:}{Right side of interval;}
#'   \item{\code{cens}:}{Indicator whether interval resulted from censoring or truncation
#'   (1 = censoring, 0 = truncation);}
#'   \item{\code{id}:}{(required) Identifier for the observation this interval belongs to (numeric/integer);}
#' } Note that the truncation intervals need to be in the form (N, Inf] with N a numeric value.
#' @param supportdf A data from containing 2 columns indicating the support intervals:
#' \describe{
#'   \item{\code{L}:}{Left side of interval;}
#'   \item{\code{R}:}{Right side of interval;}
#' }
#' 
#' @keywords internal
#' 
#' 

existenceNPMLE <- function(intervals, supportdf){
  #Note that supportdf will come supplied already ordered according to L (and thus R).
  
  
  #-------------------------Data checks----------------------------
  #Check that all columns are present
  if(!all(c("L", "R", "cens", "id") %in% colnames(intervals))){
    stop("Please provide a data.frame with at least columns named L, R, cens and id.")
  }
  intervals <- intervals[,c("L", "R", "cens", "id")]
  
  #Check that each id has 2 observations: 1 censoring and 1 truncation
  if(!all(table(intervals[, c("cens", "id")] == 1))){
    stop("Please provide a data.frame where each observation has a corresponding censoring and truncation interval.")
  }
  
  
  #-----------------Function Body-----------------------------------
  #Give the intervals a temporary name so we can match them later
  intervals["name"] <- 1:nrow(intervals)
  intervals <- intervals[, c("name", "L", "R", "cens", "id")]
  
  #Get indices of truncation and censoring intervals
  cens_idx <- which(intervals[, "cens"] == 1)
  trunc_idx <- setdiff(1:nrow(intervals), cens_idx)
  
  #Initiate Astar and Bstar as in Hudgens (2005)
  #In these lists we will gather the indices of support intervals 
  Astar <- vector(mode = "list", length = length(cens_idx))
  Bstar <- vector(mode = "list", length = length(trunc_idx))
  
  #Determine Astar
  cens_intervals <- intervals[cens_idx, c("L", "R", "id")]
  #An interval [L_1, R_1] \subseteq [L_2, R_2] if:
  #L_1 >= L2 and R_1 <= R_2
  RsmallerR <- outer(supportdf[,"R"], cens_intervals$R, "<=")
  LgreaterL <- outer(supportdf[,"L"], cens_intervals$L, ">=")
  #Determine intersections (both have to be true, so just sum, if sum==2 then intersect)
  intersec <- (RsmallerR + LgreaterL)
  #We get the row and column indices to construct the graph
  for(i in 1:ncol(intersec)){
    Astar[[i]] <- which(intersec[,i] == 2)
  }

  
  #Determine Bstar
  trunc_intervals <- intervals[trunc_idx, c("L", "R", "id")]
  #An interval [L_1, R_1] \subseteq [L_2, Inf] if:
  #L_1 >= L_2
  LgreaterL <- outer(supportdf[,"L"], trunc_intervals$L, ">=")
  for(i in 1:ncol(LgreaterL)){
    Bstar[[i]] <- which(LgreaterL[,i])
  }
  
  #Create graph to check existence
  rowcol_idx <- data.frame()
  Anames <- intervals[cens_idx, "id"]
  Bnames <- intervals[trunc_idx, "id"]
  for(i in 1:length(Astar)){
    for(j in 1:length(Bstar)){
      if(i != j){
        if(all(Astar[[i]] %in% Bstar[[j]])){
          rowcol_idx <- rbind(rowcol_idx, c(Anames[i], Bnames[j]))
        }
      }
    }
  }
  if(nrow(rowcol_idx) == 0){
    g = make_empty_graph(n = 1, directed = TRUE)
    mle_exist = TRUE
  } else{
    g <- graph_from_data_frame(rowcol_idx, directed = TRUE, vertices = cens_intervals[, c("id", "L", "R")])
    #NPMLE exists if components are strongly connected
    mle_exist <- is_connected(g, mode = "strong")
  }
  
  return(list(dir_graph = g,
              exist_mle = mle_exist))
}









