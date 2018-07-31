#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Bayesian Variable Selection
//' 
//' Produces a draw of coefficients with a Minnesota Prior.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param Z Matrix containing the time series of the explanatory variables.
//' @param A Matirx.
//' @param Gamma sfg.
//' @param Sigma_i The inverse of the current draw of the variance-covariance matrix.
//' @param pos_res fg.
//' @param lpr_prior_0 Numeric vector containing the prior mean of the coefficients.
//' @param lpr_prior_1 Covariance matrix of the coefficients.
//' 
//' @return Vector of parameter values.
//' 
//' @references
//' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection. \emph{Journal of Applied Econometrics}, 28(2), 204--230.
//' 
// [[Rcpp::export]]
arma::mat bvs(arma::mat y, arma::mat Z, arma::mat A, arma::mat Gamma, arma::mat Omega_i, arma::vec pos_res, arma::vec lpr_prior_0, arma::vec lpr_prior_1)
{
  int n = y.n_rows;
  int t = y.n_cols;
  int nomega = Omega_i.n_rows;
  
  bool const_par = true;
  if (A.n_cols > 1){
    const_par = false;
  }
  
  bool sv = false;
  if (nomega > n) {
    sv = true;
  }
  
  arma::mat S_i = arma::zeros<arma::mat>(n * t, n * t);
  if (sv) {
    for (int i = 0; i < t; i++) {
      S_i.submat(i * n, i * n, (i + 1) * n - 1, (i + 1) * n - 1 ) = Omega_i.rows(i * n, (i + 1) * n - 1);
    }
  } else {
    S_i = kron(Omega_i, arma::eye<arma::mat>(t, t));  
  }
  arma::vec yvec = vectorise(y);
  arma::vec g = arma::zeros<arma::vec>(1);
  arma::mat AG = Gamma * A;
  arma::mat theta0 = AG;
  arma::mat theta1 = AG;
  arma::vec l0_res = arma::zeros<arma::vec>(n * t);
  arma::vec l1_res = arma::zeros<arma::vec>(n * t);
  
  double l0 = 0;
  double l1 = 0;
  double bayes = 0;
  double bayes_rand = 0;
  
  pos_res = shuffle(pos_res) - 1;
  int pos_size = size(pos_res)(0);
  int k = 0;

  for (int j = 0; j < pos_size; j++){
    k = pos_res(j);
    g = log(arma::randu<arma::vec>(1));
    if (Gamma(k, k) == 1 && g(0) > lpr_prior_0(k)){continue;}
    if (Gamma(k, k) == 0 && g(0) > lpr_prior_1(k)){continue;}
    if ((Gamma(k, k) == 1 && g(0) < lpr_prior_0(k)) || (Gamma(k, k) == 0 && g(0) < lpr_prior_1(k))){
      theta0 = AG;
      theta1 = AG;
      theta0.row(k) = A.row(k);
      if (const_par) {
        theta1.row(k) = 0;
        l0_res = yvec - Z * theta0;
        l1_res = yvec - Z * theta1;
      } else {
        //t0 = vectorise(theta0);
        //t1 = vectorise(theta1);
        theta1.row(k) = arma::zeros<arma::mat>(1, t);
        for (int i = 0; i < t; i++){
          l0_res.rows(i * n, (i + 1) * n - 1) = yvec.rows(i * n, (i + 1) * n - 1) - Z.rows(i * n, (i + 1) * n - 1) * theta0.col(i);
          l1_res.rows(i * n, (i + 1) * n - 1) = yvec.rows(i * n, (i + 1) * n - 1) - Z.rows(i * n, (i + 1) * n - 1) * theta1.col(i);
        }
      }
      l0 = -.5 * arma::as_scalar(trans(l0_res) * S_i * l0_res) + arma::as_scalar(lpr_prior_0(k));
      l1 = -.5 * arma::as_scalar(trans(l1_res) * S_i * l1_res) + arma::as_scalar(lpr_prior_1(k));
      bayes = l0 - l1;
      bayes_rand = log(arma::as_scalar(arma::randu<arma::vec>(1)));
      if (bayes >= bayes_rand){
        Gamma(k, k) = 1;
      } else {
        Gamma(k, k) = 0;
      }
    }
  }
  return Gamma;
}