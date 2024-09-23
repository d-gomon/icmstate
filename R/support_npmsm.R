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

