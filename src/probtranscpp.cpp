#include <iostream>
#include <Rcpp.h>
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <string>
#include <unordered_map>
using namespace Rcpp;
using namespace Eigen;



//' Extract intensity matrices from msfit object
//' 
//' @param object List containing msfit object with elements: Haz (DataFrame) and trans (IntegerMatrix)
//' @return List containing intensity_matrices (3D array) and unique_times
// [[Rcpp::export]]
List get_intensity_matrices_cpp(const List& object) {
  
  List haz_data = object["Haz"];
  IntegerMatrix trans_matrix = object["trans"];
  
  NumericVector time_vec = haz_data["time"];
  IntegerVector trans_vec = haz_data["trans"];
  NumericVector haz_vec = haz_data["Haz"];
  
  // Get unique times
  std::vector<double> unique_times_vec(time_vec.begin(), time_vec.end());
  std::sort(unique_times_vec.begin(), unique_times_vec.end());
  unique_times_vec.erase(std::unique(unique_times_vec.begin(), unique_times_vec.end()), 
                         unique_times_vec.end());
  
  int n_times = unique_times_vec.size();
  int n_states = trans_matrix.nrow();
  
  // Create transition mapping (from -> to)
  std::vector<std::pair<int, int>> trans_map;
  for(int i = 0; i < n_states; i++) {
    for(int j = 0; j < n_states; j++) {
      if(trans_matrix(i, j) > 0) {
        trans_map.push_back({i, j});
      }
    }
  }
  
  int n_trans = trans_map.size();
  
  // Initialize intensity matrices
  NumericVector intensity_matrices(n_states * n_states * n_times);
  std::fill(intensity_matrices.begin(), intensity_matrices.end(), 0.0);
  
  // Process each transition
  for(int k = 0; k < n_trans; k++) {
    int from = trans_map[k].first;
    int to = trans_map[k].second;
    int trans_id = k + 1; // R uses 1-based indexing
    
    // Get hazard data for this transition
    std::vector<double> haz_for_trans;
    for(int i = 0; i < haz_vec.size(); i++) {
      if(trans_vec[i] == trans_id) {
        haz_for_trans.push_back(haz_vec[i]);
      }
    }
    
    // Calculate intensities as differences
    std::vector<double> intensities(n_times);
    intensities[0] = haz_for_trans[0];
    for(int i = 1; i < n_times; i++) {
      intensities[i] = haz_for_trans[i] - haz_for_trans[i-1];
    }
    
    // Assign to 3D array (stored as 1D with manual indexing)
    for(int t = 0; t < n_times; t++) {
      int idx = from + to * n_states + t * n_states * n_states;
      intensity_matrices[idx] = intensities[t];
    }
  }
  
  // Set diagonal elements (1 - row sums)
  for(int t = 0; t < n_times; t++) {
    for(int i = 0; i < n_states; i++) {
      double row_sum = 0.0;
      for(int j = 0; j < n_states; j++) {
        if(i != j) {
          int idx = i + j * n_states + t * n_states * n_states;
          row_sum += intensity_matrices[idx];
        }
      }
      int diag_idx = i + i * n_states + t * n_states * n_states;
      intensity_matrices[diag_idx] = 1.0 - row_sum;
    }
  }
  
  // Set dimensions
  intensity_matrices.attr("dim") = IntegerVector::create(n_states, n_states, n_times);
  
  return List::create(
    Named("intensity_matrices") = intensity_matrices,
    Named("unique_times") = NumericVector(unique_times_vec.begin(), unique_times_vec.end())
  );
}





//' Calculate transition probabilities (main function)
//' 
//' @param int_mat List from get_intensity_matrices_cpp
//' @param predt Prediction time
//' @param direction String: "forward" or "fixedhorizon"
//' @param as_df Boolean: return as dataframe format
//' @return Array or List of transition probabilities
//' @import RcppEigen
// [[Rcpp::export]]
SEXP probtrans_D_cpp(List int_mat, double predt, std::string direction = "forward", bool as_df = false) {
 
 NumericVector intensity_matrices = int_mat["intensity_matrices"];
 NumericVector unique_times = int_mat["unique_times"];
 IntegerVector dims = intensity_matrices.attr("dim");
 
 int n_states = dims[0];
 int n_times = dims[2];
 
 // Find relevant matrices and times
 std::vector<int> relevant_idx;
 std::vector<double> times;
 
 if(direction == "forward") {
   for(int i = 0; i < n_times; i++) {
     if(predt < unique_times[i]) {
       relevant_idx.push_back(i);
     }
   }
   times.push_back(predt);
   for(int idx : relevant_idx) {
     times.push_back(unique_times[idx]);
   }
 } else { // fixedhorizon
   for(int i = 0; i < n_times; i++) {
     if(predt >= unique_times[i]) {
       relevant_idx.push_back(i);
     }
   }
   for(int idx : relevant_idx) {
     times.push_back(unique_times[idx]);
   }
 }
 
 int n_matrices = relevant_idx.size();
 if(n_matrices == 0) {
   stop("Prediction time 'predt' must be within study time (boundaries not allowed)");
 }
 
 // Initialize probability matrix P as identity
 MatrixXd P = MatrixXd::Identity(n_states, n_states);
 
 // Calculate output dimensions
 int n_out_times;
 if(direction == "forward") {
   n_out_times = n_matrices + 1;
 } else {
   n_out_times = (std::find(unique_times.begin(), unique_times.end(), predt) != unique_times.end()) ? 
   n_matrices + 1 : n_matrices + 2;
 }
 
 // Initialize output array
 NumericVector out(n_states * (n_states + 1) * n_out_times);
 
 if(direction == "forward") {
   // Set times in first column
   for(int j = 0; j < n_states; j++) {
     for(int t = 0; t < n_out_times; t++) {
       int idx = j + 0 * n_states + t * n_states * (n_states + 1);
       out[idx] = times[t];
     }
   }
   
   // Set initial probabilities (identity matrix)
   for(int i = 0; i < n_states; i++) {
     for(int j = 0; j < n_states; j++) {
       int idx = i + (j + 1) * n_states + 0 * n_states * (n_states + 1);
       out[idx] = P(i, j);
     }
   }
   
   // Forward multiplication
   for(int t = 0; t < n_matrices; t++) {
     MatrixXd Q = MatrixXd::Zero(n_states, n_states);
     for(int i = 0; i < n_states; i++) {
       for(int j = 0; j < n_states; j++) {
         int mat_idx = i + j * n_states + relevant_idx[t] * n_states * n_states;
         Q(i, j) = intensity_matrices[mat_idx];
       }
     }
     P = P * Q;
     
     // Store result
     for(int i = 0; i < n_states; i++) {
       for(int j = 0; j < n_states; j++) {
         int idx = i + (j + 1) * n_states + (t + 1) * n_states * (n_states + 1);
         out[idx] = P(i, j);
       }
     }
   }
 } else { // fixedhorizon
   // Set times
   bool predt_in_times = std::find(unique_times.begin(), unique_times.end(), predt) != unique_times.end();
   
   for(int j = 0; j < n_states; j++) {
     if(predt_in_times) {
       for(int t = 0; t < n_out_times; t++) {
         int idx = j + 0 * n_states + t * n_states * (n_states + 1);
         out[idx] = (t == 0) ? 0.0 : times[t - 1];
       }
     } else {
       for(int t = 0; t < n_out_times; t++) {
         int idx = j + 0 * n_states + t * n_states * (n_states + 1);
         if(t == 0) out[idx] = 0.0;
         else if(t <= n_matrices) out[idx] = times[t - 1];
         else out[idx] = predt;
       }
     }
   }
   
   // Set final probabilities
   int final_idx = predt_in_times ? n_matrices : n_matrices + 1;
   for(int i = 0; i < n_states; i++) {
     for(int j = 0; j < n_states; j++) {
       int idx = i + (j + 1) * n_states + final_idx * n_states * (n_states + 1);
       out[idx] = P(i, j);
     }
   }
   
   if(!predt_in_times) {
     for(int i = 0; i < n_states; i++) {
       for(int j = 0; j < n_states; j++) {
         int idx = i + (j + 1) * n_states + (n_matrices + 1) * n_states * (n_states + 1);
         out[idx] = P(i, j);
       }
     }
   }
   
   // Backward multiplication
   for(int t = 0; t < n_matrices; t++) {
     MatrixXd Q = MatrixXd::Zero(n_states, n_states);
     int mat_idx_pos = n_matrices - 1 - t;
     for(int i = 0; i < n_states; i++) {
       for(int j = 0; j < n_states; j++) {
         int mat_idx = i + j * n_states + relevant_idx[mat_idx_pos] * n_states * n_states;
         Q(i, j) = intensity_matrices[mat_idx];
       }
     }
     P = Q * P;
     
     // Store result
     int out_idx_pos = n_matrices - t;
     for(int i = 0; i < n_states; i++) {
       for(int j = 0; j < n_states; j++) {
         int idx = i + (j + 1) * n_states + out_idx_pos * n_states * (n_states + 1);
         out[idx] = P(i, j);
       }
     }
   }
 }
 
 // Set dimensions
 out.attr("dim") = IntegerVector::create(n_states, n_states + 1, n_out_times);
 
 if(as_df) {
   // Convert to list of data frames (similar to R output)
   List res(n_states);
   for(int s = 0; s < n_states; s++) {
     NumericMatrix df_mat(n_out_times, n_states + 1);
     for(int t = 0; t < n_out_times; t++) {
       for(int c = 0; c < n_states + 1; c++) {
         int idx = s + c * n_states + t * n_states * (n_states + 1);
         df_mat(t, c) = out[idx];
       }
     }
     CharacterVector col_names(n_states + 1);
     col_names[0] = "time";
     for(int i = 1; i <= n_states; i++) {
       col_names[i] = "pstate" + std::to_string(i);
     }
     DataFrame df = DataFrame(df_mat);
     df.names() = col_names;
     res[s] = df;
   }
   return res;
 } else {
   // Return permuted array
   NumericVector res(n_out_times * (n_states + 1) * n_states);
   for(int t = 0; t < n_out_times; t++) {
     for(int c = 0; c < n_states + 1; c++) {
       for(int s = 0; s < n_states; s++) {
         int old_idx = s + c * n_states + t * n_states * (n_states + 1);
         int new_idx = t + c * n_out_times + s * n_out_times * (n_states + 1);
         res[new_idx] = out[old_idx];
       }
     }
   }
   res.attr("dim") = IntegerVector::create(n_out_times, n_states + 1, n_states);
   return res;
 }
}
 
