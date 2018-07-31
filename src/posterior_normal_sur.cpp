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