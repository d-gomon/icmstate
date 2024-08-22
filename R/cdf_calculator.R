#' Determine the cdf of the transitions in an illness-death model in the 
#' Frydman(1995) setting.
#' 
#' @description Only used in \code{\link{msm_frydman}} to determine support sets.
#' 
#' @keywords internal
#' @noRd


cdf_calculator <- function(z_lambda, supportMSM, data_idx){
  
  #Estimator of F12(x)
  cumsum_z_I <- cumsum(z_lambda$z[1:supportMSM$I])
  F12 <- function(x){
    #Find index of final right bound that is smaller equal than value
    right_bound_idx <- Position(isTRUE, x >= supportMSM$Q_mat[,2], right = TRUE)
    if(x < supportMSM$Q_mat[1,1]){
      #If evaluation value smaller than smallest left bound, return 0
      return(0)
    } else if(is.na(right_bound_idx)){
      #If smaller than first upper bound, but larger than first lower bound
      return(cumsum_z_I[1])
    } else if(right_bound_idx != supportMSM$I){
      #if any right bound is larger than value:
      if(x < supportMSM$Q_mat[right_bound_idx+1, 1]){
        #if next left bound is larger than value
        return(cumsum_z_I[right_bound_idx])
      } else{
        #else value is between the next interval in Q_i
        return(cumsum_z_I[c(right_bound_idx, right_bound_idx +1)])
      }
    } else{
      #Else we are beyond the largest bound, so return maximal value
      return(cumsum_z_I[right_bound_idx])
    }
  }
  out <- list(F12 = F12)
  
  #Estimator of F13(x)
  if(data_idx$K != 0){
    cumsum_z_I_dash <- cumsum(z_lambda$z[(supportMSM$I+1):(supportMSM$I + data_idx$K)])
    F13 <- function(x){
      #Find index of failure time that is smaller equal than value
      right_bound_idx <- Position(isTRUE, x >= data_idx$e_k_star, right = TRUE)
      if(x < data_idx$e_k_star[1]){
        return(0)
      } else{
        return(cumsum_z_I_dash[right_bound_idx])  
      }
    }
    out <- c(out, F13 = F13)
  }
  
  #Estimator of Lambda23(x)
  cumsum_lambda <- unname(cumsum(z_lambda$lambda))
  Lambda23 <- function(x){
    right_bound_idx <- Position(isTRUE, x >= data_idx$t_n_star, right = TRUE)
    if(x < data_idx$t_n_star[1]){
      return(0)
    } else{
      return(cumsum_lambda[right_bound_idx])
    }
  }
  out <- c(out, Lambda23 = Lambda23)
  
  F23 <- function(x){
    warning("This function is not working correctly at the moment: returning exp(-H(t)) instead of Prodi 1-dH(t). Should be fixed later.")
    return(1 - exp(-Lambda23(x)))
  }
  out <- c(out, F23 = F23)
  
  return(out)
}