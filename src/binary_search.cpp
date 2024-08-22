#include <Rcpp.h>
using namespace Rcpp;


//' Binary search - Larger or equal to
//' 
//' @description Find index of first value larger or equal to x in a sorted array
//' 
//' @keywords internal
//' @param v - sorted vector instance
//' @param data - value to search
//' @return index of first value larger or equal than data, -1 otherwise
// [[Rcpp::export]]
int binary_search_larger_equal(NumericVector v, double data) {
  auto it = std::lower_bound(v.begin(), v.end(), data);
  if (it == v.end()) {
    return -1;
  }  else {
    std::size_t index = std::distance(v.begin(), it);
    return index;
  }
}


//' Binary search - Larger 
//' 
//' @description Find index of first value larger than x in a sorted array
//' 
//' @keywords internal
//' @param v - sorted vector instance
//' @param data - value to search
//' @return index of first value larger than data -1, -1 otherwise
// [[Rcpp::export]]
int binary_search_larger(NumericVector v, double data) {
  auto it = std::upper_bound(v.begin(), v.end(), data);
  if (it == v.end()) {
    return -1;
  } else {
    std::size_t index = std::distance(v.begin(), it);
    return index;
  }
}

