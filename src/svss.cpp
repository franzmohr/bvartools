#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Stochastic Search Variable Selection
//' 
//' \code{ssvs} employs stochastic search variable selection as proposed by George et al. (2008)
//' to produce a draw of the precision matrix of the coefficients in a VAR model.
//' 
//' @param a an M-dimensional vector of coefficient draws.
//' @param tau0 an M-dimensional vector of prior standard deviations for restricted
//' coefficients in vector \code{a}.
//' @param tau1 an M-dimensional vector of prior standard deviations for unrestricted
//' coefficients in vector \code{a}.
//' @param prob_prior an M-dimensional vector of prior inclusion probabilites for the coefficients
//' in vector \code{a}.
//' @param include an integer vector specifying the positions of coefficients in vector \code{a}, which should be
//' included in the SSVS algorithm. If \code{NULL} (default), SSVS will be applied to all coefficients.
//' 
//' @details The function employs stochastic search variable selection (SSVS) as proposed
//' by George et al. (2008) to produce a draw of the diagonal inverse prior covariance matrix
//' \eqn{\underline{V}^{-1}} and the corresponding vector of inclusion parameters \eqn{\lambda}
//' of the vectorised coefficient matrix \eqn{a = vec(A)} for the VAR model
//' \deqn{y_t = A x_t + u_t,}
//' where \eqn{y_{t}} is a K-dimensional vector of endogenous variables,
//' \eqn{x_{t}} is a vector of explanatory variabes
//' and the error term is \eqn{u_t \sim \Sigma}.
//' 
//' @return A named list containing two components:
//' \item{v_i}{an \eqn{M \times M} inverse prior covariance matrix.}
//' \item{lambda}{an M-dimensional vector of inclusion parameters.}
//' 
//' @examples
//' # Prepare data
//' data("e1")
//' data <- diff(log(e1))
//' temp <- gen_var(data, p = 2, deterministic = "const")
//' y <- temp$Y
//' x <- temp$Z
//' k <- nrow(y)
//' tt <- ncol(y)
//' m <- k * nrow(x)
//' 
//' # OLS estimates for semiautomatic approach
//' ols <- tcrossprod(y, x) %*% solve(tcrossprod(x))
//' # OLS error covariance matrix
//' sigma_ols <- tcrossprod(y - ols %*% x) / (tt - nrow(x))
//' # Covariance matrix
//' cov_ols <- kronecker(solve(tcrossprod(x)), sigma_ols)
//' # Standard errors
//' se_ols <- matrix(sqrt(diag(cov_ols)))
//' 
//' # SSVS priors
//' tau0 <- se_ols * 0.1 # If excluded
//' tau1 <- se_ols * 10 # If included
//' 
//' # Prior for inclusion parameter
//' prob_prior <- matrix(0.5, m)
//' 
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(c(tau1^2), m)
//' 
//' # Initial value of Sigma
//' sigma_i <- solve(sigma_ols)
//' 
//' # Draw parameters
//' a <- post_normal(y = y, x = x, sigma_i = sigma_i,
//'                  a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//' 
//' # Run SSVS
//' lambda <- ssvs(a = a, tau0 = tau0, tau1 = tau1,
//'                prob_prior = prob_prior)
//' 
//' @references
//' 
//' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
//' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
//' \url{https://doi.org/10.1016/j.jeconom.2007.08.017}
//' 
// [[Rcpp::export]]
  Rcpp::List ssvs(arma::vec a, arma::vec tau0, arma::vec tau1, arma::vec prob_prior,
                Rcpp::Nullable<Rcpp::IntegerVector> include = R_NilValue) {
  int m = a.n_rows;
  arma::vec tau0_sq = arma::square(tau0);
  arma::vec tau1_sq = arma::square(tau1);
  arma::vec u0 = 1 / tau0 % arma::exp(-(arma::square(a) / (2 * tau0_sq))) % (1 - prob_prior);
  arma::vec u1 = 1 / tau1 % arma::exp(-(arma::square(a) / (2 * tau1_sq))) % prob_prior;
  arma::vec p_post = u1 / (u0 + u1);
  arma::vec lambda = arma::ones<arma::vec>(m);
  arma::mat V = arma::zeros<arma::mat>(m, m);
  V.diag() = 1 / tau1_sq;
  double temp;
  
  arma::vec pos_res(m);
  for (int l = 0; l < m; l++) {
    pos_res(l) = l;
  }
  arma::uvec ex;
  if (include.isNotNull()) {
    ex = Rcpp::as<arma::uvec>(include) - 1;
    pos_res = pos_res.elem(ex);
  }
  
  pos_res = shuffle(pos_res);
  int pos_size = size(pos_res)(0);
  int k;
  for (int i = 0; i < pos_size; i++){
    k = pos_res(i);
    temp = Rcpp::as<double>(Rcpp::rbinom(1, 1, p_post(k)));
    lambda(k) = temp;
    if (temp == 0) {
      V(k, k) = 1 / tau0_sq(k);
    }
  }
  return Rcpp::List::create(Rcpp::Named("v_i") = V,
                            Rcpp::Named("lambda") = lambda);
}