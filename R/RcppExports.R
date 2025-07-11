# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Binary search - Larger or equal to
#' 
#' @description Find index of first value larger or equal to x in a sorted array
#' 
#' @keywords internal
#' @param v - sorted vector instance
#' @param data - value to search
#' @return index of first value larger or equal than data, -1 otherwise
binary_search_larger_equal <- function(v, data) {
    .Call(`_icmstate_binary_search_larger_equal`, v, data)
}

#' Binary search - Larger 
#' 
#' @description Find index of first value larger than x in a sorted array
#' 
#' @keywords internal
#' @param v - sorted vector instance
#' @param data - value to search
#' @return index of first value larger than data -1, -1 otherwise
binary_search_larger <- function(v, data) {
    .Call(`_icmstate_binary_search_larger`, v, data)
}

#' Extract intensity matrices from msfit object
#' 
#' @param object List containing msfit object with elements: Haz (DataFrame) and trans (IntegerMatrix)
#' @return List containing intensity_matrices (3D array) and unique_times
get_intensity_matrices_cpp <- function(object) {
    .Call(`_icmstate_get_intensity_matrices_cpp`, object)
}

#' Calculate transition probabilities (main function)
#' 
#' @param int_mat List from get_intensity_matrices_cpp
#' @param predt Prediction time
#' @param direction String: "forward" or "fixedhorizon"
#' @param as_df Boolean: return as dataframe format
#' @return Array or List of transition probabilities
#' @import RcppEigen
probtrans_D_cpp <- function(int_mat, predt, direction = "forward", as_df = FALSE) {
    .Call(`_icmstate_probtrans_D_cpp`, int_mat, predt, direction, as_df)
}

