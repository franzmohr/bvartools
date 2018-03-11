#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Bayesian Variable Selection
//' 
//' Produces a draw of selection variables using the pseudo-althorithm from Korobilis (2013).
//' 
//' @param y a $n x T$ matrix containing the time series of the dependent variable.
//' @param Z Matrix containing the time series of the explanatory variables.
//' @param A Matirx.
//' @param Gamma sfg.
//' @param Sigma_i The inverse of the current draw of the variance-covariance matrix.
//' @param pos_res fg.
//' @param lpr_prior_0 Numeric vector containing the prior mean of the coefficients.
//' @param lpr_prior_1 Covariance matrix of the coefficients.
//' 
//' @return a vector of inclusion indicators.
//' 
//' @references
//' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection. \emph{Journal of Applied Econometrics}, 28(2), 204--230.
//' 
// [[Rcpp::export]]
arma::vec bvs(arma::mat y, arma::mat Z, arma::mat A, arma::mat Gamma, arma::mat Sigma_i, arma::vec pos_res, arma::vec lpr_prior_0, arma::vec lpr_prior_1)
{
  int n = y.n_rows;
  int t = y.n_cols;
  int nsig = Sigma_i.n_rows;

  bool const_par = true;
  if (A.n_cols > 1){
    const_par = false;
  }
  
  bool sv = false;
  if (nsig > n) {
    sv = true;
  }
  
  arma::vec g = arma::zeros<arma::vec>(1);
  
  arma::mat AG = diagmat(Gamma) * A;
  arma::mat theta0 = AG;
  arma::mat theta1 = AG;
  arma::vec l0_res = arma::zeros<arma::vec>(n);
  arma::vec l1_res = arma::zeros<arma::vec>(n);
  
  double L0 = 0;
  double L1 = 0;
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
    if (Gamma(k, 0) == 1 && g(0) > lpr_prior_0(k)){continue;}
    if (Gamma(k, 0) == 0 && g(0) > lpr_prior_1(k)){continue;}
    if ((Gamma(k, 0) == 1 && g(0) < lpr_prior_0(k)) || (Gamma(k, 0) == 0 && g(0) < lpr_prior_1(k))){
      theta0 = AG;
      theta1 = AG;
      theta0.row(k) = A.row(k);
      theta1.row(k) = 0;
      L0 = 0;
      L1 = 0;
      for (int i = 0; i < t; i++){
        if (const_par){
          l0_res = y.col(i) - Z.rows(i * n, (i + 1) * n - 1) * theta0;
          l1_res = y.col(i) - Z.rows(i * n, (i + 1) * n - 1) * theta1;
        } else {
          theta0(k, i) = A(k, i + 1);
          theta1(k, i) = 0;
          l0_res = y.col(i) - Z.rows(i * n, (i + 1) * n - 1) * theta0.col(i);
          l1_res = y.col(i) - Z.rows(i * n, (i + 1) * n - 1) * theta1.col(i);
        }
        if (sv){
          L0 = L0 + as_scalar(trans(l0_res) * Sigma_i.rows(i * n, (i + 1) * n - 1) * l0_res);
          L1 = L1 + as_scalar(trans(l1_res) * Sigma_i.rows(i * n, (i + 1) * n - 1) * l1_res);
        } else {
          L0 = L0 + as_scalar(trans(l0_res) * Sigma_i * l0_res);
          L1 = L1 + as_scalar(trans(l1_res) * Sigma_i * l1_res);
        }
      }
    }
    l0 = -.5 * L0 + arma::as_scalar(lpr_prior_0(k));
    l1 = -.5 * L1 + arma::as_scalar(lpr_prior_1(k));
    bayes = l0 - l1;
    bayes_rand = log(as_scalar(arma::randu<arma::vec>(1)));
    if (bayes >= bayes_rand){
      Gamma(k, 0) = 1;
    } else {
      Gamma(k, 0) = 0;
    }
  }
  return Gamma; 
}