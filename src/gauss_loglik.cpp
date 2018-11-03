#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Log-Likelihood
//' 
//' Calculates the log likelihood of the residuals of a multivariate time series.
//' 
//' @param y an \eqn{n x T} matrix of residuals.
//' @param Sigma an \eqn{n x n} or \eqn{nT x n} variance-covariance matrix.
//' @param Sigma_i the inverse of Sigma.
//' 
// [[Rcpp::export]]
arma::vec gauss_loglik(arma::mat y, arma::mat Sigma, arma::mat Sigma_i) {
  int n = y.n_rows;
  int t = y.n_cols;
  double d = 0;
  double tpi = pow(2 * arma::datum::pi, -.5 * n);
  
  arma::vec result = arma::zeros<arma::vec>(t);
  
  if (Sigma_i.n_rows > n) {
    for (int i = 0; i < t; i++){
      d = pow(det(Sigma.rows(i * n, (i + 1) * n - 1)), -.5);
      result(i) = tpi * d * exp(-.5 * as_scalar(trans(y.col(i)) * Sigma_i.rows(i * n, (i + 1) * n - 1) * y.col(i)));
    }
  } else {
    d = pow(det(Sigma), -.5);
    for (int i = 0; i < t; i++){
      result(i) = tpi * d * exp(-.5 * as_scalar(trans(y.col(i)) * Sigma_i * y.col(i)));
    }
  }
  return result;
}