#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw from Normal Distribution
//' 
//' Produces a draw from a posterior density of a Gaussian process and normal prior.
//' 
//' @param y \eqn{n x T} matrix containing the time series of the dependent variable.
//' @param x \eqn{k x T} matrix containing the time series of the explanatory variables.
//' @param Sigma_i inverse of the \eqn{n x n} variance-covariance matrix.
//' @param bprior \eqn{m x 1} numeric vector containing the prior mean of the coefficients.
//' @param Vprior_i inverse of the \eqn{m x m} covariance matrix of the coefficients.
//' 
//' @details The function produces a vectorised draw of the \eqn{n x k} coefficient matrix \eqn{A} of the model
//' 
//' \deqn{y_{t} = A x_{t} a + \epsilson_{t},}
//' 
//' where \eqn{y_{t}} is a \eqn{n x 1} vector of endogenous variables in period \eqn{t}, \eqn{x_{t}} is a \eqn{k x 1}
//' vector of explanatory variables, and the error term \eqn{\epsilon_{t}} is normally distributed with zero mean
//' and variance-covariance matrix \eqn{\Sigma}. \eqn{k} is the number of explanatory variables and \eqn{m} it the
//' total number of estimated coefficients.
//' 
//' @return a vector of coefficient draws.
//' 
// [[Rcpp::export]]
arma::vec posterior_normal(arma::mat y, arma::mat x, arma::mat Sigma_i, arma::vec bprior, arma::mat Vprior_i) {
  int nvars = y.n_rows * x.n_rows;
  
  arma::mat Vpost = arma::inv(Vprior_i + arma::kron(x * arma::trans(x), Sigma_i));
  arma::vec bpost = Vpost * (Vprior_i * bprior + vectorise(Sigma_i * y * arma::trans(x)));

  arma::mat U;
  arma::mat V;
  arma::vec s;
  arma::svd(U, s, V, Vpost);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  
  arma::vec z = arma::randn<arma::vec>(nvars);
  arma::vec result = bpost + A * z;
  return result;
}