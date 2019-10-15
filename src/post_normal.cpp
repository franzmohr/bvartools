#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw from a Normal Distribution
//' 
//' Produces a draw of coefficients from a normal posterior density.
//' 
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param x an \eqn{M \times T} matrix of explanatory variables.
//' @param sigma_i the inverse of the \eqn{K \times K} variance-covariance matrix.
//' @param a_prior a \eqn{KM \times 1} numeric vector of prior means.
//' @param v_i_prior the inverse of the \eqn{KM \times KM} prior covariance matrix.
//' 
//' @details The function produces a vectorised posterior draw \eqn{a} of the
//' \eqn{K \times M} coefficient matrix \eqn{A} for the model
//' \deqn{y_{t} = A x_{t} + u_{t},}
//' where \eqn{y_{t}} is a K-dimensional vector of endogenous variables,
//' \eqn{x_{t}} is an M-dimensional vector of explanatory variabes
//' and the error term is \eqn{u_t \sim \Sigma}.
//' 
//' For a given prior mean vector \eqn{\underline{a}} and prior covariance matrix \eqn{\underline{V}}
//' the posterior covariance matrix is obtained by
//' \deqn{\overline{V} = \left[ \underline{V}^{-1} + \left(X X^{\prime} \otimes \Sigma^{-1} \right) \right]^{-1}}
//' and the posterior mean by
//' \deqn{\overline{a} = \overline{V} \left[ \underline{V}^{-1} \underline{a} + vec(\Sigma^{-1} Y X^{\prime}) \right],}
//' where \eqn{Y} is a \eqn{K \times T} matrix of the endogenous variables and \eqn{X} is an \eqn{M \times T} matrix of
//' the explanatory variables.
//' 
//' @examples
//' # Prepare data
//' data("e1")
//' data <- diff(log(e1))
//' temp <- gen_var(data, p = 2, deterministic = "const")
//' y <- temp$Y
//' x <- temp$Z
//' k <- nrow(y)
//' t <- ncol(y)
//' m <- k * nrow(x)
//' 
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(0.1, m)
//' 
//' # Initial value of inverse Sigma
//' sigma_i <- solve(tcrossprod(y) / t)
//' 
//' # Draw parameters
//' a <- post_normal(y = y, x = x, sigma_i = sigma_i,
//'                  a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//' 
//' @return A vector.
//' 
//' @references
//' 
//' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
//' 
// [[Rcpp::export]]
arma::vec post_normal(arma::mat y, arma::mat x, arma::mat sigma_i, arma::vec a_prior, arma::mat v_i_prior) {
  int m = y.n_rows * x.n_rows;
  int n_mu = a_prior.n_rows;
  int n_v = v_i_prior.n_rows;
  if (m != n_mu) {
    Rcpp::stop("Argument 'a_prior' does not contain the required amount of elements.");
  }
  if (m != n_v) {
    Rcpp::stop("Argument 'v_i_prior' does not contain the required amount of elements.");
  }
  
  arma::mat v_post = arma::inv(v_i_prior + arma::kron(x * arma::trans(x), sigma_i));
  arma::vec a_post = v_post * (v_i_prior * a_prior + vectorise(sigma_i * y * arma::trans(x)));

  arma::mat U;
  arma::vec s;
  arma::eig_sym(s, U, v_post);

  return a_post + (U * arma::diagmat(sqrt(s)) * arma::trans(U)) * arma::randn<arma::vec>(m);
}