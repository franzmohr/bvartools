// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export(.draw_forecast)]]
arma::mat draw_forecast(int &i, // index of draw
                        int &k, // number of endogenous vars
                        int &p, // number of lags of endogenous vars
                        arma::mat &a0_i, // A0
                        bool &use_a,
                        arma::mat &a_, // A
                        arma::mat &sigma, // Sigma
                        arma::mat pred) { // Data matrix for prediction
  
  const int n_ahead = pred.n_cols - 1;
  int n_tot = 0;
  arma::mat a;
  if (use_a) {
    n_tot = a_.n_cols / k;
    a = arma::reshape(a_.row(i - 1), k, n_tot);
  }
  arma::vec eigval, u;
  arma::mat eigvec;
  
  // Forecast iterations
  for (int j = 0; j < n_ahead; j++) {
    
    // Generate random error
    arma::eig_sym(eigval, eigvec, arma::reshape(sigma.row(i - 1), k, k));
    u = eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t() * arma::randn(k);
    
    // Forecast for next period
    if (use_a) {
      if (p == 0) {
        pred.submat(0, j + 1, k - 1, j + 1) = a0_i * a * pred.submat(k, j, k + n_tot - 1, j) + a0_i * u; 
      } else {
        pred.submat(0, j + 1, k - 1, j + 1) = a0_i * a * pred.col(j) + a0_i * u;  
      }
    } else {
      pred.submat(0, j + 1, k - 1, j + 1) = pred.col(j) + a0_i * u; 
    }
    
    // Update lags of endogenous variables
    if (p > 1) {
      for (int l = 0; l < (p - 1); l++) {
        pred.submat((l + 1) * k, j + 1, (l + 2) * k - 1, j + 1) = pred.submat(l * k, j, (l + 1) * k - 1, j); 
      }
    }
  }
  
  return pred;
}