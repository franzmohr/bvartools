#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Seemingly Unrelated Regression
//' 
//' Produces a draw of normally distributed coefficients with using SUR.
//' 
//' @param y a matrix containing the time series of the dependent variable.
//' @param Sigma_i a \eqn{n \times n} or \eqn{n * T \times n} (diagonal) variance matrix. 
//' @param A0_mu_prior a numeric vector containing the prior mean of the coefficients.
//' @param A0_V_i_prior a matrix containing the prior precisions of the coefficients.
//' @param tvp Matrix containing the time series of the explanatory variables.
//' 
//' @return Vector of parameter values.
//' 
// [[Rcpp::export]]
Rcpp::List draw_structural_parameters(arma::mat y, arma::mat Sigma_i, arma::mat A0_mu_prior,
                                      arma::mat A0_V_i_prior, bool tvp) {
  
  int n = y.n_rows;
  int t = y.n_cols;
  int p1 = 0;
  int p2 = 0;
  bool sv = Sigma_i.n_rows > n;
  
  arma::mat y_temp;
  arma::mat x_temp;
  arma::mat Sigma_i_temp;
  arma::mat A0_mu_prior_temp;
  arma::mat A0_V_i_prior_temp;
  arma::mat A0_mu_post;
  arma::mat A0_V_post;
  arma::mat U;
  arma::mat V;
  arma::vec s;
  arma::mat A0 = arma::eye<arma::mat>(n, n);
  arma::mat A0_temp;
  arma::mat y_dot = y * 0;
  arma::mat S_i;
  arma::mat ZHi;
  arma::mat ZHZ;
  arma::mat ZHy;
  
  if (tvp) {
    
  } else {
    for (int i = 1; i < n; i++) {
      p1 = i * (i - 1) / 2;
      p2 = i * (i + 1) / 2 - 1;
      y_temp = y.row(i);
      x_temp = -y.rows(0, i - 1);
      A0_mu_prior_temp = A0_mu_prior.rows(p1, p2);
      A0_V_i_prior_temp = A0_V_i_prior.submat(p1, p1, p2, p2);
      if (sv) {
        ZHi = arma::zeros<arma::mat>(i, t);
        ZHZ = arma::zeros<arma::mat>(i, i);
        ZHy = arma::zeros<arma::mat>(i, 1);
        S_i = arma::zeros<arma::mat>(t, t);
        for (int j = 0; j < t; j++){
          S_i.diag()(j) = Sigma_i(j * n + i, i);
        }
        //ZHi = arma::trans(x_temp) * S_i;
        ZHi = x_temp * S_i;
        //ZHZ = ZHi * x_temp;
        ZHZ = ZHi * arma::trans(x_temp);
        //ZHy = ZHi * arma::vectorise(y_temp);
        ZHy = ZHi * arma::trans(y_temp);
        A0_V_post = arma::inv(A0_V_i_prior_temp + ZHZ);
        A0_mu_post = A0_V_post * (A0_V_i_prior_temp * A0_mu_prior_temp + ZHy);
      } else {
        Sigma_i_temp = Sigma_i(i,i);
        A0_V_post = arma::inv(A0_V_i_prior_temp + arma::kron(x_temp * arma::trans(x_temp), Sigma_i_temp));
        A0_mu_post = A0_V_post * (A0_V_i_prior_temp * A0_mu_prior_temp + vectorise(Sigma_i_temp * y_temp * arma::trans(x_temp)));
      }
      arma::svd(U, s, V, A0_V_post);
      A0_temp = A0_mu_post + U * arma::diagmat(sqrt(s)) * arma::trans(V) * arma::randn<arma::vec>(i);
      A0.submat(i, 0, i, i - 1) = arma::trans(A0_temp); 
    }
    y_dot = A0 * y;
  }
  
  return Rcpp::List::create(A0, y_dot);
}