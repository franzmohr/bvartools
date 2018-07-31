#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Draw Cointegration Matrix with Minnesota Prior
//' 
//' Produces a draw of coefficients of the cointegration matrix with a Minnesota Prior.
//' 
//' @param y Matrix containing the time series of the dependent variable.
//' @param ect Matrix containing the time series of the explanatory variables.
//' @param alpha Loading matrix of Pi.
//' @param Sigma_i Inverse covariance matrix.
//' 
//' @return List of the matrix \eqn{\Pi} and the cointegration matrix \eqn{\beta}.
//' 
//' @references
//' Koop, G., Leon-Gonzalez, R., & Strachan, R. W. (2010). Efficient posterior simulation for cointegrated models with priors on the cointegration space. \emph{Econometric Reviews}, 29(2), 224--242.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List posterior_cointegration(arma::mat y, arma::mat ect, arma::mat beta, arma::mat Sigma_i, arma::mat G_i, double v_i, arma::mat P_i) {
  int n = y.n_rows;
  int nect = ect.n_rows;
  int r = beta.n_cols;
  int nalpha = n * r;
  int nbeta = nect * r;

  arma::mat x = arma::trans(beta) * ect;
  arma::mat Alpha = arma::zeros<arma::mat>(n, r);
  
  arma::mat U;
  arma::vec s;
  arma::mat V;

  arma::mat alpha_V_post = arma::inv(arma::kron(v_i * arma::trans(beta) * P_i * beta, G_i) + arma::kron(x * arma::trans(x), Sigma_i));
  arma::vec alpha_mu_post = alpha_V_post * arma::vectorise(Sigma_i * y * arma::trans(x));
  arma::svd(U, s, V, alpha_V_post);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  
  arma::vec z = arma::randn<arma::vec>(nalpha);
  arma::mat alpha = reshape(alpha_mu_post + A * z, n, r);
  
  arma::svd(U, s, V, alpha);
  Alpha = U.cols(0, r - 1) * arma::trans(V);

  arma::mat Beta_V_post = arma::inv(arma::kron(arma::trans(Alpha) * G_i * Alpha, v_i * P_i) + arma::kron(arma::trans(Alpha) * Sigma_i * Alpha, ect * arma::trans(ect)));
  arma::vec Beta_mu_post = Beta_V_post * vectorise(ect * arma::trans(y) * Sigma_i * Alpha);
  
  arma::svd(U, s, V, Beta_V_post);
  A = U * arma::diagmat(sqrt(s)) * arma::trans(V);
  z = arma::randn<arma::vec>(nbeta);
  arma::mat Beta = arma::reshape(Beta_mu_post + A * z, nect, r);

  arma::svd(U, s, V, Beta);
  alpha = Alpha * V * arma::diagmat(s) * arma::trans(V);
  beta = U.cols(0, r - 1) * arma::trans(V);
  
  arma::mat Pi = alpha * arma::trans(beta);
  Rcpp::List result = Rcpp::List::create(Pi, alpha, beta);
  return result;
}