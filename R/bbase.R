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
  n <- dim(P)[2]
  D <- diff_D(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B
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

#' Non-lagged Differences (adjusted)
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


diff_D <- function(x, differences = 1L, ...){
  xlen <- dim(x)[1L]
  if (1L * differences >= xlen)
    return(x[0L]) # empty, but of proper mode
  i1 <- -1L
  for (i in seq_len(differences))
    x <- x[i1, , drop = FALSE] -
  x[-nrow(x):-(nrow(x)), , drop = FALSE]
  x
}