#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draw from Wishart Distribution
//' 
//' Produces a posterior draw from a Wishart Distribution.
//' 
//' @param y a \eqn{n x T} matrix of the times series.
//' @param df_post posterior degrees of freedom.
//' @param V_prior \eqn{n x n} prior variance-covariance matrix.
//' 
//' @return a matrix.
//' 
// [[Rcpp::export]]
arma::mat post_wishart(arma::mat y, int df_post, arma::mat V_prior) {
  arma::mat V_i_post = arma::inv(V_prior + y * arma::trans(y));
  arma::mat result = arma::wishrnd(V_i_post, df_post);
  return result;
}