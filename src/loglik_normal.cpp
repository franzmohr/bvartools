#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculates the log-likelihood of a multivariate normal distribution.
//' 
//' @param u a \eqn{K \times T} matrix of residuals.
//' @param sigma a \eqn{K \times K} or \eqn{KT \times K} variance-covariance matrix.
//' 
//' @details The log-likelihood is calculated for each vector in period \eqn{t} as
//' \deqn{-\frac{K}{2} \ln 2\pi - \frac{1}{2} \ln |\Sigma_t| -\frac{1}{2} u_t^\prime \Sigma_t^{-1} u_t},
//' where \eqn{u_t = y_t - \mu_t}.
//' 
//' @examples
//' 
//' # Load data
//' data("e1")
//' e1 <- diff(log(e1))
//' 
//' # Generate VAR model
//' data <- gen_var(e1, p = 2, deterministic = "const")
//' y <- t(data$data$Y)
//' x <- t(data$data$Z)
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
  int k = u.n_rows;
  int t = u.n_cols;
  
  arma::mat dmat = arma::eye(k, k);
  arma::vec result = arma::zeros<arma::vec>(t);
  double part_b, part_c;
  double part_a = -k * log(2 * arma::datum::pi) / 2 ;
  if (sigma.n_rows > u.n_rows) {
    for (int i = 0; i < t; i++){
      part_b = -log(arma::det(sigma.rows(i * k, (i + 1) * k - 1))) / 2;
      part_c = -arma::as_scalar(trans(u.col(i)) * arma::solve(sigma.rows(i * k, (i + 1) * k - 1), dmat) * u.col(i)) / 2;
      result(i) = part_a + part_b + part_c;
    }
  } else {
    part_b = -log(arma::det(sigma)) / 2;
    for (int i = 0; i < t; i++){
      part_c = -arma::as_scalar(trans(u.col(i)) * arma::solve(sigma, dmat) * u.col(i)) / 2;
      result(i) = part_a + part_b + part_c;
    }
  }
  return result;
}