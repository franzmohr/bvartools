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
// [[Rcpp::export]]
Rcpp::List posterior_pi(arma::mat y, arma::mat ect , arma::mat alpha, arma::mat Sigma_i) {
  int nect = ect.n_rows;
  int r = alpha.n_cols;
  int nvars = nect*r;
  arma::mat Alpha = alpha*sqrtmat_sympd(inv(trans(alpha)*alpha));
  arma::mat Vpost = inv(kron(trans(Alpha)*Sigma_i*Alpha,ect*trans(ect)));
  arma::vec bpost = Vpost * vectorise(ect*trans(y)*Sigma_i*Alpha);
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval,eigvec,Vpost);
  arma::mat A = eigvec * diagmat(sqrt(eigval));
  arma::vec z = arma::randn<arma::vec>(nvars);
  arma::mat Beta = reshape(bpost + A * z, nect, r);
  arma::mat kappa = sqrtmat_sympd(trans(Beta)*Beta);
  arma::mat beta = Beta*inv(kappa);
  alpha = Alpha*kappa;
  arma::mat Pi = alpha*trans(beta);
  Rcpp::List result = Rcpp::List::create(Pi, beta);
  return result; 
}