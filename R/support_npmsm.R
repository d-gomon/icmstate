#' Numerically find the support of the transitions from a converged npmsm algorithm
#' 
#' @description For each transition in \code{tmat}, determine all consecutive 
#' bins with non-zero (higher than \code{cutoff}) transition intensities.
#' These then determine the numerical support of the transition.
#' 
#' 
#' @param npmsm Output from \code{\link{npmsm}} function or an \code{\link[mstate:msfit]{msfit}} object.
#' @param cutoff Above which value is a mass in a bin considered to be non-zero? 
#' Default = 1e-8. Note that this is independent of bin size, so can be tricky!!
#' 
#' @importFrom mstate to.trans2 
#'
#' @export
#' 
#' @returns A list containing a list for each transition. Each transition 
#' specific list contains the support intervals for that transition in a matrix 
#' with 3 named columns L, R and dA, indicating the left/right-endpoints of the 
#' support intervals and the change in the estimated intensities over this 
#' support interval. 
#' 
#' @examples 
#' require(mstate)
#' require(ggplot2)
#' #Generate from an illness-death model with exponential transitions with 
#' #rates 1/2, 1/10 and 1 for 10 subjects over a time grid.
#' gd <- sim_weibmsm(tmat = trans.illdeath(), shape = c(1,1,1),
#'                   scale = c(2, 10, 1), n_subj = 10, obs_pars = c(2, 0.5, 20), 
#'                   startprobs = c(0.9, 0.1, 0))
#' #Fit 2 models: 1 with at most 4 iterations and 1 with at most 20
#' mod1 <- npmsm(gd, trans.illdeath(), maxit = 4)
#' 
#' #Determine support numerically:
#' mod1_supp <- support_npmsm(mod1)
#' mod1_supp[[1]]




support_npmsm <- function(npmsm, cutoff = 1e-8){
  
  #Remove CRAN notes
  trans <- NULL
  
  #Calculate number of transitions
  if(inherits(npmsm, "npmsm")){
    tmat <- npmsm$tmat
    tmat2 <- npmsm$tmat2
    A <- npmsm$A
  } else if(inherits(npmsm, "msfit")){
    tmat <- npmsm$trans
    tmat2 <- to.trans2(tmat)
    A <- npmsm
  }
  
  H <- nrow(tmat) # no of states
  M <- nrow(tmat2) # no of transitions
  
  supportHein <- vector(mode = "list", length = M)
  for(j in 1:M){
    #Subset on transition
    temptrans <- subset(A$Haz, trans == j)
    pointmasses <- c(temptrans$Haz[1], diff(temptrans$Haz))
    
    #Determine support for current transition
    supportHeini <- matrix(ncol = 3, nrow = 0)
    Ltemp <- NULL
    Rtemp <- NULL
    dAtemp <- 0
    for(i in 1:length(pointmasses)){
      dAi <- pointmasses[i] #Store the current "pointmass" A(t_{i}) - A(t_{i-1})
      if(dAi > cutoff && is.null(Ltemp)){
        if(i == 1){
          Ltemp <- 0
        } else{
          Ltemp <- temptrans$time[i-1]  
        }
      }
      if(!is.null(Ltemp)){
        if(dAi > cutoff){
          dAtemp <- dAtemp + dAi  
        }
      }
      if(i == length(pointmasses) && dAi > cutoff){
        Rtemp <- temptrans$time[i]
        supportHeini <- rbind(supportHeini, c(Ltemp, Rtemp, dAtemp))
        Ltemp <- NULL
        dAtemp <- 0
      }
      if(dAi <= cutoff && !is.null(Ltemp)){
        Rtemp <- temptrans$time[i-1]
        supportHeini <- rbind(supportHeini, c(Ltemp, Rtemp, dAtemp))
        Ltemp <- NULL
        dAtemp <- 0
      }
    }
    colnames(supportHeini) <- c("L", "R", "dA")
    supportHein[[j]]$support <- supportHeini
    supportHein[[j]]$transition <- unlist(c(tmat2[match(j, tmat2$transno), c("from", "to")]))
  }
  names(supportHein) <- tmat2$transname
  return(supportHein)
  
}

