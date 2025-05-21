#' Compute a B-spline basis 
#' 
#' Copied from icpack/bases.R, modified diff() function for speed.
#'
#' @export
#' @param x The vector of values for which the basis is to be evaluated
#' @param xl The left boundary of the domain
#' @param xr The right boundary of the domain
#' @param nseg The number of inter-knot segments on the domain
#' @param bdeg The degree of the B-splines (2 means quadratic, 3 means cubic, and so on)
#' @return A matrix containing the basis
#'
#'
#' @keywords internal
#' @noRd
#' 
#' @examples
#' x = runif(100)
#' B = bbase_D(x, 0, 1, 20, 3)

bbase_D <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) {
  # Construct B-spline basis
  dx <- (xr - xl) / nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower_D, bdeg)
  n <- length(knots)
  D <- diff_D(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B
}


#' Compute a B-spline basis for a single time point
#' 
#' Similar to bbase_D, but sped up for single time points.
#'
#' @export
#' @param x The value for which the basis is to be evaluated (ONLY SINGLE VALUE ALLOWED)
#' @param xl The left boundary of the domain
#' @param xr The right boundary of the domain
#' @param nseg The number of inter-knot segments on the domain
#' @param bdeg The degree of the B-splines (2 means quadratic, 3 means cubic, and so on)
#' @return A vector containing the basis
#'
#'
#' @keywords internal
#' @noRd
#' 
#' @examples
#' x = 0.02
#' B = bbase_singletime(x, 0, 1, 20, 3)




bbase_singletime <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3){
  #Future speed-up ideas:
  #Pre-calculate knots outside this function. Then index only over the indices of interest
  #Pre-calculate D = diff_D(diag(n), differences = bdeg + 1). 
  #Essentially re-calculating the same matrix every time, especially as n is basically the same every time
  
  #Only works if x is a single time - not longer for vectors!
  
  
  #Compare compute speed to old implementation
  #microbenchmark(new = bbase_singletime(0.09, 0, 1, nseg = 20, bdeg = 3),
  #               old = bbase_D(0.09, 0, 1, nseg = 20, bdeg = 3),
  #               times = 500)
  
  #Knot distances
  dx <- (xr - xl)/nseg
  
  #Knot locations: we only need to calculate B-spline values within the segment we are interested in
  #For this (see Appendix C.1) we only need the left and right side of segment x is located in
  #and 3 extra knots to each side to calculate the tpower at the knots
  #Left side of segment x is located in: xl + floor((x-xl)/dx)
  #Right side: xl + ceiling((x-xl)/dx)
  x_segment <- (x-xl)/dx
  floor_x_segment <- floor(x_segment)
  ceil_x_segment <- ceiling(x_segment)
  
  start <- floor_x_segment - bdeg
  end <- ceil_x_segment + bdeg
  
  #Calculate knots
  #knots <- seq(xl + (floor_x_segment - bdeg) * dx, 
  #             xl + (ceil_x_segment + bdeg) * dx, 
  #             by = dx)
  knots <- xl + dx * seq.int(start, end)
  #Faster than using seq()
  
  #Left and Right side of segment which x is located in + bdeg knots to both sides (for tpower calculations)
  #P <- outer(x, knots, tpower_D, bdeg)
  P <- (x - knots)^bdeg * (x >= knots)
  #Can directly apply tpower_D function as x is single value.
  #Now we need to scale the contributions
  #n <- length(knots)
  n <- end - start + 1L 
  #No need to re-calculate D every time. This function is run extremely often
  t_D <- diff_D_t(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  #D <- diff_D(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t_D
  
  #Need to determine the indices which we now have to fill in.
  #We start from the right side of the segment
  #We know that only n_knots-bdeg-1 entries are ever non-zero
  
  B_out <- vector(mode = "numeric", length = nseg + bdeg) #The output should have a column for every spline.
  length_fill <- n - bdeg -1 #Number of elements to fill
  fill_indices <- ceiling(x_segment) + seq_len(length_fill) - 1 #seq_len is faster than 0:length_fill
  if(floor_x_segment == ceil_x_segment){
    fill_indices <- fill_indices + 1
  }
  B_out[fill_indices] <- B
  B_out
}


bbase_singletime_cached <- function(x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3, fix_pars){
  #Same as bbase_singletime, but the following quantities are cached (pre-calculated):
  #- dx, t_D, 
  #Note that there are 2 distinct t_D matrices
  #If floor_x_segment == ceil_x_segment, t_D has one dimension less!
  #We call this t_D_trunc
  #This is either calculated using diff_D_t(diag(n)) with n = 2*bdeg + 1 or 2*bdeg + 2
  
  #Knot distances
  dx <- fix_pars[["bbase_dx"]]
  
  #knot locations
  x_segment <- (x-xl)/dx
  floor_x_segment <- floor(x_segment)
  ceil_x_segment <- ceiling(x_segment)
  
  start <- floor_x_segment - bdeg
  end <- ceil_x_segment + bdeg
  
  knots <- xl + dx * seq.int(start, end)
  #Faster than using seq()
  
  #Apply tpower_D to knots
  P <- (x - knots)^bdeg * (x >= knots)
  #n <- length(knots)
  n <- end - start + 1L 
  #No need to re-calculate D every time. This function is run extremely often
  if (floor_x_segment != ceil_x_segment) {
    t_D <- fix_pars[["bbase_t_D"]]
  } else {
    t_D <- fix_pars[["bbase_t_D_trunc"]]
  }
  #D <- diff_D(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- fix_pars[["bbase_const"]] * P %*% t_D
  
  #Fill output vector
  B_out <- vector(mode = "numeric", length = nseg + bdeg) #The output should have a column for every spline.
  length_fill <- n - bdeg -1 #Number of elements to fill
  fill_indices <- ceiling(x_segment) + seq_len(length_fill) - 1 #seq_len is faster than 0:length_fill
  if(floor_x_segment == ceil_x_segment){
    fill_indices <- fill_indices + 1
  }
  B_out[fill_indices] <- B
  B_out
}



#' Truncated power function
#' 
#' Same as JOPS::tpower(), internally defined for quick reference
#' 
#' @keywords internal
#' @noRd
#' 
#' 

tpower_D <- function(x, t, p) {
  # Truncated p-th power function
  (x - t) ^ p * (x >= t)
}

#' Non-lagged (transposed) Differences (adjusted)
#' 
#' Adjusted version of base diff() function, can only be used with matrices and 
#' only with lag = 1L
#' 
#' @param x A matrix containing the values to be differenced
#' @param differences An integer indicating the order of the difference
#' @param ... Further arguments to diff
#' 
#' @return A matrix with the difference operation carried out on each column separately
#' 
#' 
#' @keywords internal
#' @noRd
#' 
#' 


diff_D <- function(x, differences = 1L) {
  for (i in seq_len(differences)) {
    x <- x[-1L, , drop = FALSE] - x[-nrow(x), , drop = FALSE]
  }
  x
}

#' @keywords internal

diff_D_t <- function(x, differences = 1L){
  for(i in seq_len(differences)){
    x <- x[, -1L, drop = FALSE] - x[, -ncol(x), drop = FALSE]
  }
  x
}


