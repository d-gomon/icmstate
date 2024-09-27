#' @title icmstate
#' @description Non-parametric estimation of transition intensities in 
#' interval-censored multi-state models
#' 
#' 
#' @details
#' Allows for the estimation of transition intensities in interval-censored
#' multi-state models using the approach of Gomon and Putter (2024) (Multinomial likelihood) or 
#' Gu et al. (2023) (Poisson likelihood).
#'
#' @author Maintainer: Daniel Gomon <dgstatsoft@gmail.com> [aut, cre] \cr
#' Hein Putter [aut]
#' 
#' @name icmstate-package
#' @aliases icmstate-package
#' @docType package
#' @references Y. Gu, D. Zeng, G. Heiss, and D. Y. Lin, 
#' Maximum likelihood estimation for semiparametric regression models with 
#' interval-censored multistate data, Biometrika, Nov. 2023, doi:10.1093/biomet/asad073
#' 
#' D. Gomon and H. Putter, 
#' Non-parametric estimation of transition intensities in interval censored Markov multi-state models without loops,
#' arXiv, 2024, doi:10.48550/arXiv.2409.07176
#' 
#' 
#' 
#' 
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib icmstate, .registration = TRUE
## usethis namespace: end
"_PACKAGE"