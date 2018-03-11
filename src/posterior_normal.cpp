#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//' Posterior Draw
//' 
//' Produces a draw from a posterior density of a Gaussion process with normal prior.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param x Matrix containing the time series of the explanatory variables.
//' @param Sigma_i Inverse matrix of the covariance matrix of a Gaussian process.
//' @param bprior Numeric vector with prior mean of the coefficients.
//' @param Vprior_i Inverse prior covariance matrix of the coefficients.
//' 
//' @details
//' \code{y} and \code{x} must have the form \eqn{n x T} and \eqn{m x T}, respectively,
//' where \eqn{n} is the number of endogenous variables, \eqn{m} the number of explanatory variables and
//' \eqn{T} the number of observations.
//' 
//' @return Vector of parameter values.
//' 
// [[Rcpp::export]]
arma::vec posterior_normal(arma::mat y, arma::mat x, arma::mat Sigma_i, arma::vec bprior, arma::mat Vprior_i) {
  int nvars = y.n_rows * x.n_rows;
  arma::mat Vpost = inv(Vprior_i + kron(x*trans(x),Sigma_i));
  arma::vec bpost = Vpost*(Vprior_i*bprior + vectorise(Sigma_i*y*trans(x)));

  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval,eigvec,Vpost);
  arma::mat A = eigvec * diagmat(sqrt(eigval));
  arma::vec z = arma::randn<arma::vec>(nvars);
  arma::vec result = bpost + A * z;
  return result;
}