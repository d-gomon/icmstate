#' @title Check if closed interval is contained in other closed interval
#'
#' @description Function which takes as input two matrices with 2 columns each 
#' and checks whether each interval in the first matrix is contained within 
#' each interval in the second matrix.
#'
#' @param A Two column matrix containing intervals to be checked for being 
#' contained in B
#' @param B Two column matrix containing intervals possibly overlapping 
#' the intervals in A
#' @param B.left.open Are the intervals in B left-open?
#' @param B.right.open Are the intervals in B right-open?
#'
#'
#' @return Matrix of size (nrow(A) * nrow(B)) with binary values indicating 
#' whether the intervals in A are contained in B
#'




AsubsetB <- function(A, B, B.left.open = FALSE, B.right.open = FALSE){
  #Store final information
  nrow_A <- nrow(A)
  nrow_B <- nrow(B)
  
  if(isTRUE(B.left.open)){
    `%left.operator%` <- function(e1, e2) e1 > e2
  } else{
    `%left.operator%` <- function(e1, e2) e1 >= e2
  }
  
  if(isTRUE(B.right.open)){
    `%right.operator%` <- function(e1, e2) e1 < e2
  } else{
    `%right.operator%` <- function(e1, e2) e1 <= e2
  }
  out <- matrix(NA, nrow = nrow_A, ncol = nrow_B)
  
  for(i in 1:nrow_A){
    for(j in 1:nrow_B){
      if(A[i,1] %left.operator% B[j,1] & A[i,2] %right.operator% B[j,2]){
        out[i, j] <- 1
      } else{
        out[i,j] <- 0
      } 
    }
  }
  rownames(out) <- rownames(A)
  colnames(out) <- rownames(B)
  out
}


#' @title Check if closed interval is contained in half-open infinite interval
#'
#' @description Function which takes as input a matrix with 2 columns and a 
#' vector indicating left points of intervals [b, Infinity)
#' and checks whether each interval in the matrix is contained within 
#' the corresponding interval derived from b.
#'
#' @param A Two column matrix containing intervals to be checked for being 
#' contained in b
#' @param b Vector indicating left point in corresponding [b, Infinity) interval
#'
#'
#' @return Matrix of size (nrow(A) * length(b)) with binary values indicating 
#' whether the intervals in A are contained in the ones induced by b
#'


Alargerb <- function(A, b){
  #Store final information
  nrow_A <- nrow(A)
  length_B <- length(b)
  
  out <- matrix(NA, nrow = nrow_A, ncol = length_B)
  
  for(i in 1:nrow_A){
    for(j in 1:length_B){
      if(A[i,1] > b[j]){
        out[i, j] <- 1
      } else{
        out[i,j] <- 0
      } 
    }
  }
  rownames(out) <- rownames(A)
  colnames(out) <- names(b)
  out
}


#' @title Check if event time is contained within half-open interval
#'
#' @description Function which takes as input a vector a with event times 
#' and a 2 column matrix B representing half-open intervals (l, R] and checks whether 
#' each event time is contained in each half-open interval.
#'
#' @param a Vector of event times
#' @param B Two column matrix containing intervals
#'
#'
#' @return Matrix of size (length(a) * nrow(B)) with binary values indicating 
#' whether the event times in a are contained in the intervals of B
#'

ainB <- function(a, B){
  length_a <- length(a)
  nrow_B <- nrow(B)
  
  out <- matrix(NA, nrow = length_a, ncol = nrow_B)
  
  for(i in 1:length_a){
    for(j in 1:nrow(B)){
      if(a[i] > B[j,1] & a[i] <= B[j,2]){
        out[i, j] <- 1
      } else{
        out[i,j] <- 0
      } 
    }
  }
  rownames(out) <- names(a)
  colnames(out) <- rownames(B)
  out
}


#' @title Check if event time is larger/equal than other event time
#'
#' @description Function which takes as input a vector a with event times 
#' and a vector b with event times and checks whether each entry in a is 
#' larger or equal than the entries in b.
#'
#' @param a Vector of event times
#' @param b Vector of event times
#'
#'
#' @return Matrix of size (length(a) * length(b)) with binary values indicating 
#' whether the event times in a are larger than the event times in b
#'
ageqb <- function(a, b){
  length_a <- length(a)
  length_b <- length(b)
  
  out <- matrix(NA, nrow = length_a, ncol = length_b)
  
  for(i in 1:length_a){
    for(j in 1:length_b){
      if(a[i] >= b[j]){
        out[i, j] <- 1
      } else{
        out[i,j] <- 0
      } 
    }
  }
  rownames(out) <- names(a)
  colnames(out) <- names(b)
  out
}

#' @title Check if event time is larger than other event time
#'
#' @description Function which takes as input a vector a with event times 
#' and a vector b with event times and checks whether each entry in a is 
#' larger  than the entries in b.
#'
#' @param a Vector of event times
#' @param b Vector of event times
#'
#'
#' @return Matrix of size (length(a) * length(b)) with binary values indicating 
#' whether the event times in a are larger than the event times in b
#'
agreaterb <- function(a, b){
  length_a <- length(a)
  length_b <- length(b)
  
  out <- matrix(NA, nrow = length_a, ncol = length_b)
  
  for(i in 1:length_a){
    for(j in 1:length_b){
      if(a[i] > b[j]){
        out[i, j] <- 1
      } else{
        out[i,j] <- 0
      } 
    }
  }
  rownames(out) <- names(a)
  colnames(out) <- names(b)
  out
}



#' @title Check if half-open intervals intersect with event times
#'
#' @description Function which takes as input a 2 column matrix of half-open 
#' intervals and a vector of event times and returns the event times that 
#' intersect with the half-open intervals.
#'
#' @param A Two column matrix containing intervals
#' @param b Vector of event times
#' @param A.left.open Are the intervals in A open on the left side? Default = FALSE.
#'
#'
#' @return Numeric vector of event times from b that intersect with A.
#'


Aintersectb <- function(A, b, A.left.open = FALSE){
  nrow_A <- nrow(A)
  length_b <- length(b)
  
  if(isTRUE(A.left.open)){
    `%left.operator%` <- function(e1, e2) e1 > e2
  } else{
    `%left.operator%` <- function(e1, e2) e1 >= e2
  }
  
  out <- numeric(0)
  
  for(i in 1:length_b){
    for(j in 1:nrow_A){
      if(b[i] %left.operator% A[j,1] & b[i] <= A[j,2]){
        out <- c(out, b[i])
        break
      }
    }
  }
  out
}



#' @title Calculate the product of intensities over interval decided by 
#' failure times
#' 
#' @description This function calculates \eqn{\prod_{\lambda}}{prod(\lambda)} G as defined in 
#' Frydman (1995), with \eqn{t_n}{t_n} the failure times. Note that length(t_n) must be 
#' equal to length(lambda)
#'
#' @param lambda Intensities of the 2->3 transition
#' @param t_n Unique failure times, same length as lambda
#' @param Q Matrix (2 column) containing support intervals as rows
#' @param A Matrix (2 column) containing censoring intervals as rows
#'
#'
#'


prod_lambda_G_base <- function(lambda, t_n, Q, A){
  #Function to calculate Prod_\lambda {G}
  #lambda is input vector of length N (current estimate of lambda)
  #z is input vector of length I (current estimate of z)
  #G is a half-open interval (r_i, R_m] (vector of length 2)
  #t_n_star vector of length N associated with lambda
  #N = length(lambda) = length(t_n_star)
  N <- length(t_n)
  out <- 1
  I <- nrow(Q)
  M <- nrow(A)
  
  for(j in 1:N){
    if(t_n[j] > A[1] && t_n[j] <= A[2]){
      out <- out * (1-lambda[j])
    }
  }
  return(out)
}





#' Find submatrices of an array that are the identity matrix
#'
#' @description Given an array of stacked block matrices, find in the direction of the \code{dim}'th dimension
#' all the matrices that are the identity matrix.
#'
#' @param arr An (dim_block, dim_block, r) dimensional array, with r any positive 
#' integer. The dimensions can be any permutation of above display.
#' @param dim The dimension in which to perform the search (of size r above)
#' @param dim_block The size of the block dimension
#'
#'
#'
#' @keywords internal
#' @noRd
#'



which_identity_arr <- function(arr, dim, dim_block){
  apply(arr, dim, function(x) !any(diag(x, names = FALSE) != 1))
}





