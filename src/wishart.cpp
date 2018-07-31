#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//' Wishart Distribution
//' 
//' Produces a single draw from the Wishart distribution.
//' 
//' @param V a positive definite scale matrix.
//' @param n an integer specifying the degrees of freedom.
//' 
//' @return a matrix.
//' 
// [[Rcpp::export]]
arma::mat wishart(arma::mat M, int n) {
  arma::mat U;
  arma::mat V;
  arma::vec s;
  svd(U, s, V, M);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(V) * arma::randn<arma::vec>(M.n_rows);
  A = A * arma::trans(A);
  return A;
}