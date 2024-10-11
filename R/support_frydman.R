#' @title Find support for Illness-death model. Only 
#' applicable to Frydman (1995) setting!!
#' 
#' @details Finds the support sets Q for the interval censored transition
#' 1 -> 2. The transitions 2->3 and 1->3 are exactly observed in this setting, 
#' therefore do not require determination of the support.
#' 
#' 
#' @param data_idx List containing all the necessary data, output from the 
#' \code{\link{msm_frydman}} function.
#' 
#' @return A: Matrix containing the intervals indicating when the transitions 
#' were censored.
#' unique_L_times: Union of failure times in Group 4 and censoring times in 
#' Group 1.
#' L_tilde: Left intervals for constructing the Q intervals
#' R_tilde: Right intervals for constructing the Q intervals
#' Q_mat: Matrix containing the support set for the 1->2 transition
#' I: number of intervals in the support set Q.
#' 
#'
#' @references Frydman, H. (1995). Nonparametric Estimation of a Markov 
#' 'Illness-Death' Process from Interval- Censored Observations, with 
#' Application to Diabetes Survival Data. Biometrika, 82(4), 773-789. 
#' \doi{10.2307/2337344} 
#'
#'
#' @keywords internal
#'


support_frydman <- function(data_idx){
  #Given processed MSM data from msm_frydman(), calculates the support of the 
  #cdf. 
  
  #-----------------Some data preparation----------------------#
  #Unique (failure times in Group 4) union (censoring times in Group 1)
  unique_L_times <- unique(c(data_idx$matdata_g4[, 5], data_idx$matdata_g1[, 5]))
  #A = all half-open intervals where 1->2 transitions have been censored
  A = data_idx$matdata_g34[, c(3,4), drop = FALSE]
  #Order A by left intervals.
  A <- A[order(A[, 1]), , drop = FALSE]
  
  #-----------------Calculate A \cap (T \cup D_{(0, 0)})----------------------#
  A_intersect_TunionD <- Aintersectb(A, unique_L_times, A.left.open = TRUE)
  
  if(isTRUE(data_idx$ltrunc)){
    #------------Determine v_r \in \overline{A}, 1 \leq r \leq N^*-------------#
    #------------Only when left truncated data------------------------#
    vr_in_overlineA <- Aintersectb(A, data_idx$v_r)
  }
  
  #--------------------Construct L_dash and R_dash-------------------------#
  #L_dash is the union of the sets on the bottom of P777 in Frydman(1995)
  L_dash <- c(unique(data_idx$matdata[, "L"]),A_intersect_TunionD)
  e_max <- ifelse(data_idx$K == 0, -1, max(data_idx$matdata_g2[, "time"]))
  R_max <- max(data_idx$matdata_g34[, "R"])
  s_max <- max(data_idx$matdata_g1[, "time"])
  #If maximal censoring time in Group 1 is larger than R_max and maximal
  #time of 1->3 transition in Group 2, add it to L_dash
  if(s_max > max(R_max, e_max)){
    L_dash <- c(L_dash, s_max)
  }
  if(isTRUE(data_idx$ltrunc)){
    R_dash <- c(unique(c(data_idx$matdata[, "R"], vr_in_overlineA)), Inf) 
  } else{
    R_dash <- c(unique(data_idx$matdata[, "R"]), Inf) 
  }
  L_dash <- sort(L_dash[!is.na(L_dash)])
  R_dash <- sort(R_dash[!is.na(R_dash)])
  
  #-------------------Calculate Q_i intervals-------------------#
  #Create matrix containing L and R_dash. 0 indicates Left bound
  #1 indicates Right bound.
  Q_mat <- matrix(c(rep(0L, length(L_dash)), rep(1L, length(R_dash)), L_dash, R_dash),
                  ncol = 2)
  #Sort Q_mat according to values
  Q_mat <- Q_mat[order(Q_mat[, 2]),]
  #Find indices where we switch from left to right bound 
  #(difference of 1 between current and next value)
  Q_mat_idx <- which(diff(Q_mat[, 1], 1) == 1)
  #Merge vectors component wise to obtain [L,R] intervals
  #Q_mat <- Q_mat[as.vector(rbind(Q_mat_idx, Q_mat_idx + 1)),]
  Q_mat <- matrix(c(Q_mat[Q_mat_idx, 2], Q_mat[Q_mat_idx + 1, 2]), ncol = 2)
  rownames(Q_mat) <- 1:nrow(Q_mat)
  colnames(Q_mat) <- c("L", "R")
  
  return(list(A = A,
              unique_L_times = unique_L_times,
              L_tilde = L_dash,
              R_tilde = R_dash,
              Q_mat = Q_mat,
              I = nrow(Q_mat)
  ))
}

