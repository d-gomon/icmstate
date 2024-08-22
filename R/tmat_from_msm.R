#'
#' Extract a \code{transMat} from a \code{msm} fit.
#' 
#'
#' @keywords internal
#' @noRd
#'



tmat_from_msm <- function(msm){
  if(!inherits(msm, "msm")){
    stop("Please provide an msm object")
  }
  
  counter <- 1
  
  imatrix <- msm$qmodel$imatrix
  for(i in 1:nrow(imatrix)){
    for(j in 1:ncol(imatrix)){
      if(imatrix[i,j] == 1){
        imatrix[i,j] <- counter
        counter <- counter + 1
      }
    }
  }
  return(imatrix)
}