#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw from Normal Distribution (SUR)
//' 
//' Produces a draw from a posterior density of a Gaussian process with normal prior using seemingly unlelated regression.
//' 
//' @param y \eqn{n x T} matrix containing the time series of the dependent variable.
//' @param Z \eqn{nT x m} matrix containing the time series of the explanatory variables.
//' @param Sigma_i inverse of the \eqn{n x n} or \eqn{n * T x n} variance-covariance matrix.
//' @param bprior \eqn{m x 1} numeric vector containing the prior mean of the coefficients.
//' @param Vprior_i inverse of the \eqn{m x m} covariance matrix of the coefficients.
//' 
//' @details The function produces a draw of the \eqn{m x 1} coefficient vector \eqn{a} of the model
//' 
//' \deqn{y_{t} = Z_{t} a + \epsilon_{t},}
//' 
//' where \eqn{y_{t}} is a \eqn{n x 1} vector of endogenous variables in period \eqn{t}, \eqn{Z_{t}} is a \eqn{n x m n}
//' matrix of explanatory variables, and the error term \eqn{\epsilon_{t}} is normally distributed with zero mean
//' and variance-covariance matrix \eqn{\Sigma_{t}}.
//' 
//' @return a vector of coefficient draws.
//' 
// [[Rcpp::export]]
arma::vec posterior_normal_sur(arma::mat y ,arma::mat Z, arma::mat Sigma_i, arma::vec bprior, arma::mat Vprior_i)
{
  int n = y.n_rows;
  int t = y.n_cols;
  int nvars = Z.n_cols;
  bool const_var = true;
  arma::mat S_i = arma::zeros<arma::mat>(n * t, n * t);
  if (Sigma_i.n_rows > n) {
    const_var = false;
  }
  
  arma::mat ZHi = arma::zeros<arma::mat>(nvars, n * t);
  arma::mat ZHZ = arma::zeros<arma::mat>(nvars, nvars);
  arma::mat ZHy = arma::zeros<arma::mat>(nvars, 1);
  
  if (const_var) {
    ZHi = arma::trans(Z) * arma::kron(arma::eye<arma::mat>(t, t), Sigma_i);
  } else {
    for (int i = 0; i < t; i++){
      S_i.submat(i * n, i * n, (i + 1) * n - 1, (i + 1) * n - 1) = Sigma_i.rows(i * n, (i + 1) * n - 1);
    }
    ZHi = arma::trans(Z) * S_i;
  }
  ZHZ = ZHi * Z;
  ZHy = ZHi * arma::vectorise(y);
  
  arma::mat Vpost = arma::inv(Vprior_i + ZHZ);
  arma::vec bpost = Vpost * (Vprior_i * bprior + ZHy);

  arma::mat U;
  arma::mat V;
  arma::vec s;
  svd(U, s, V, Vpost);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  arma::vec z = arma::randn<arma::vec>(nvars);
  arma::vec result = bpost + A * z;
  return result;
}