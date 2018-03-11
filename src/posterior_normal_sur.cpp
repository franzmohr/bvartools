#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Seemingly Unrelated Regression
//' 
//' Produces a draw of normally distributed coefficients with using SUR.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param Z Matrix containing the time series of the explanatory variables.
//' @param Sigma_i The inverse of the current draw of the variance-covariance matrix.
//' @param bprior Numeric vector containing the prior mean of the coefficients.
//' @param Vprior_i Covariance matrix of the coefficients.
//' 
//' @return Vector of parameter values.
//' 
// [[Rcpp::export]]
arma::vec posterior_normal_sur(arma::mat y ,arma::mat Z, arma::mat Sigma_i, arma::vec bprior, arma::mat Vprior_i)
{
  int n = y.n_rows;
  int t = y.n_cols;
  int nvars = Z.n_cols;
  bool sv = false;
  if (Sigma_i.n_rows>n) {
    sv = true;
  }
  arma::mat ZHi = arma::zeros<arma::mat>(nvars, n);
  arma::mat ZHZ = arma::zeros<arma::mat>(nvars, nvars);
  arma::mat ZHy = arma::zeros<arma::mat>(nvars, 1);
  if (sv) {
    for (int i=0; i < t; i++) {
      ZHi = trans(Z.rows(i * n, i * n + n - 1)) * Sigma_i.rows(i * n, i * n + n - 1);
      ZHZ = ZHZ + ZHi * Z.rows(i*n,i*n+n-1);
      ZHy = ZHy + ZHi * y.col(i);
    }
  } else {
    for (int i=0; i < t; i++) {
      ZHi = trans(Z.rows(i * n, i * n + n - 1)) * Sigma_i;
      ZHZ = ZHZ + ZHi * Z.rows(i * n, i * n + n - 1);
      ZHy = ZHy + ZHi * y.col(i);
    }
  }
  arma::mat Vpost = inv(Vprior_i + ZHZ);
  arma::vec bpost = Vpost * (Vprior_i * bprior + ZHy);

  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, Vpost);
  arma::mat A = eigvec * diagmat(sqrt(eigval));
  arma::vec z = arma::randn<arma::vec>(nvars);
  arma::vec result = bpost + A * z;
  return result;
}