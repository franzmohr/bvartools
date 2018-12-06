#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Gaussian Log-Likelihood
//' 
//' Calculates the log likelihood of a multivariate Gaussian process.
//' 
//' @param y an \eqn{n x T} matrix of residuals, where \eqn{n} is the number of
//' variables and \eqn{T} is the total amount of observations.
//' @param Sigma a constant \eqn{n x n} or time varying \eqn{nT x n} variance-covariance matrix.
//' @param Sigma_i the inverse of Sigma.
//' 
//' @return A vector of log likelihoods for each period.
//' 
// [[Rcpp::export]]
arma::vec loglik_gauss(arma::mat y, arma::mat Sigma, arma::mat Sigma_i) {
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