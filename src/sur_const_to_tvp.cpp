#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' SUR Matrix Transformation
//' 
//' Transforms a dense matrix of dimensions \eqn{KT \times M} into a sparse block
//' diagonal matrix of dimensions \eqn{KT \times MT}.
//' 
//' @param z a \eqn{KT \times M} matrix.
//' @param k integer of the number of endogenous variables.
//' @param tt integer of the number of observations.
//' 
//' @return A sparse block diagonal matrix of dimensions \eqn{KT \times MT}.
//' 
//' @examples
//' 
//' # Specify the dimensions of the dense matrix
//' k <- 2
//' tt <- 5
//' m <- 3
//' 
//' # Generate artificial data
//' z <- matrix(NA, k * tt, m)
//' for (i in 1:tt) {
//'   z[(i - 1) * k + 1:k, ] <- i
//' }
//' 
//' # Perform transformation
//' sur_const_to_tvp(z, k, tt)
//' 
//' 
// [[Rcpp::export]]
arma::sp_mat sur_const_to_tvp(const arma::mat z, const int k, const int tt) {
  
  const int m = z.n_cols;
  arma::sp_mat z_large = arma::zeros<arma::sp_mat>(k * tt, m * tt);
  for (int i = 0; i < tt; i++) {
    z_large.submat(i * k, i * m, (i + 1) * k - 1, (i + 1) * m - 1) = z.rows(i * k, (i + 1) * k - 1);
  } 
  
  return z_large;
}

/*** R

k <- 2
tt <- 5
m <- 3

z <- matrix(NA, k * tt, m)
for (i in 1:tt) {
  z[(i - 1) * k + 1:k, ] <- i
}

sur_const_to_tvp(z, k, tt)

*/