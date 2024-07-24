// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export(.draw_forecast2)]]
arma::mat draw_forecast2(int& p, // number of lags of endogenous vars
                        arma::mat& a0_i, // A0
                        bool& use_a,
                        arma::mat& a, // A
                        arma::mat& sigma, // Sigma
                        arma::mat pred) { // Data matrix for prediction
  
  int k = a0_i.n_rows;
  int n_ahead = pred.n_cols;
  int pred_end = pred.n_rows - 1;
  arma::vec eigval, u;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, sigma);
  
  // Forecast iterations
  for (int j = 0; j < n_ahead; j++) {
    
    // Generate random error
    u = eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t() * arma::randn(k);
    
    // Forecast for next period
    if (use_a) {
      pred.submat(0, j, k - 1, j) = a0_i * a * pred.submat(k, j, pred_end, j) + a0_i * u; 
    } else {
      pred.col(j) = a0_i * u; 
    }
    
    // Update lags of endogenous variables
    if (j < (n_ahead - 1) && p > 0) {
      pred.submat(k, j + 1, (p + 1) * k - 1, j + 1) = pred.submat(0, j, p * k - 1, j);  
    }
  }
  
  return pred;
}