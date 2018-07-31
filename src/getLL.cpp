#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculates the log likelihood of the residuals of a multivariate time series.
//' 
//' @param y an n x T matrix of residuals.
//' @param Sigma an n x n or nT x n variance-covariance matrix.
//' @param Sigma_i a ninverse of Sigma.
//' 
// [[Rcpp::export]]
arma::vec getLL(arma::mat y, arma::mat Sigma,arma::mat Sigma_i) {
  int n = y.n_rows;
  int t = y.n_cols;
  double d = 0;
  double tpi = pow(2 * arma::datum::pi, -.5 * n);
  //arma::vec yvec = arma::vectorise(y);
  arma::vec result = arma::zeros<arma::vec>(t);
  if (Sigma_i.n_rows > n) {
    for (int i = 0; i < t; i++){
      d = pow(det(Sigma.rows(i * n, i * n + n - 1)), -.5);
      result(i) = tpi * d * exp(-.5 * as_scalar(trans(y.col(i)) * Sigma_i.rows(i * n, i * n + n - 1) * y.col(i)));
    }
  } else {
    d = pow(det(Sigma), -.5);
    for (int i = 0; i < t; i++){
      result(i) = tpi * d * exp(-.5 * as_scalar(trans(y.col(i)) * Sigma_i * y.col(i)));
    }
  }
  return result;
}