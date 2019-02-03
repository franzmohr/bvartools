#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Stochastic Search Variable Selection
//' 
//' \code{ssvs} employs stochastic search variable selection as proposed by George et al. (2008)
//' to produce a draw of the precision matrix of the coefficients in a VAR model.
//' 
//' @param a an M-dimensional vector of coefficient draws.
//' @param tau0 an M-dimensional vector of prior standard deviations for restricted variables.
//' @param tau1 an M-dimensional vector of prior standard deviations for unrestricted variables.
//' @param prob_prior an M-dimensional vector of prior inclusion probabilites.
//' @param include an integer vector specifying the positions of variables, which should be
//' included in the SSVS algorithm. If \code{NULL} (default), SSVS will be applied to all variables.
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
//' \item{V_i}{an \eqn{M \times M} inverse prior covariance matrix.}
//' \item{lambda}{an M-dimensional vector of inclusion parameters.}
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
  return Rcpp::List::create(Rcpp::Named("V_i") = V,
                            Rcpp::Named("lambda") = lambda);
}