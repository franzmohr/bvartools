#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculates the log-likelihood of a multivariate normal distribution.
//' 
//' @param u a \eqn{K \times T} matrix of residuals.
//' @param sigma a \eqn{K \times K} or \eqn{KT \times K} variance-covariance matrix.
//' 
//' @examples
//' data("e1")
//' e1 <- diff(log(e1))
//' 
//' # Generate VAR model
//' data <- gen_var(e1, p = 2, deterministic = "const")
//' y <- data$Y
//' x <- data$Z
//'
//' # LS estimate
//' ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
//' 
//' # Residuals
//' u <- y - ols %*% x # Residuals
//'
//' # Covariance matrix
//' sigma <- tcrossprod(u) / ncol(u)
//' 
//' # Log-likelihood
//' loglik_normal(u = u, sigma = sigma)
//' 
// [[Rcpp::export]]
arma::vec loglik_normal(arma::mat u, arma::mat sigma) {
  arma::uword k = u.n_rows;
  int t = u.n_cols;
  
  arma::vec result = arma::zeros<arma::vec>(t);
  double part_b, part_c;
  double part_a = -k / 2 * log(arma::datum::pi);
  if (sigma.n_rows > k) {
    for (int i = 0; i < t; i++){
      part_b = -log(arma::det(sigma.rows(i * k, (i + 1) * k - 1))) / 2;
      part_c = -arma::as_scalar(trans(u.col(i)) * arma::inv(sigma.rows(i * k, (i + 1) * k - 1)) * u.col(i)) / 2;
      result(i) = part_a + part_b + part_c;
    }
  } else {
    part_b = -log(arma::det(sigma)) / 2;
    for (int i = 0; i < t; i++){
      part_c = -arma::as_scalar(trans(u.col(i)) * arma::inv(sigma) * u.col(i)) / 2;
      result(i) = part_a + part_b + part_c;
    }
  }
  return result;
}