#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Covariance: Vector to Matrix
//'
//' Convenience function, which takes the vector of draws of lower triangular
//' covariance coefficients and transforms it into a matrix with ones on the main
//' diagonal. In case of time varying parameters the resulting matrix will be
//' block diagonal.
//'
//' @param psi a \eqn{K (K - 1) / 2 \times 1} or \eqn{T K (K - 1) / 2 \times 1} vector of input data.
//' @param k the number \eqn{K} of endogenous variables.
//' @param tt the number \eqn{T} of observations.
//' 
//' @return A sparse, block diagonal matrix.
//' 
//' @references
//' 
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods} (2nd ed.). Cambridge: Cambridge University Press.
//' 
//' Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. \emph{The Review of Economic Studies 72}(3), 821--852. \doi{10.1111/j.1467-937X.2005.00353.x}
//' 
//' @examples
//' 
//' # Create artificial data
//' k <- 5
//' tt <- 4
//' n_covar <- (k - 1) * k / 2
//' 
//' # Constant parameters
//' psi <- matrix(1:(n_covar))
//' covar_vector_to_matrix(psi, k, tt)
//' 
//' # Time varying parameters
//' psi <- matrix(1:(n_covar * tt))
//' covar_vector_to_matrix(psi, k, tt)
//' 
// [[Rcpp::export(covar_vector_to_matrix)]]
arma::sp_mat covar_vector_to_matrix(const arma::vec psi, const int k, const int tt) {
  
  if (k < 2) {
    Rcpp::stop("Argument 'k' must be larger than 1.");
  }
  
  const arma::uword n_covar = k * (k - 1) / 2;
  arma::sp_mat psi_new;
  
  if (psi.n_elem == n_covar) {
    psi_new = arma::eye<arma::sp_mat>(k, k);
    for (int i = 1; i < k; i++) {
      psi_new.submat(i, 0, i, i - 1) = arma::trans(psi.subvec(i * (i - 1) / 2, (i + 1) * i / 2 - 1));
    }
  } else {
    psi_new = arma::eye<arma::sp_mat>(k * tt, k * tt);
    for (int j = 0; j < tt; j++) {
      for (int i = 1; i < k; i++) {
        psi_new.submat(k * j + i,
                       k * j,
                       k * j + i,
                       k * j + i - 1) = arma::trans(psi.subvec(n_covar * j + i * (i - 1) / 2,
                                                               n_covar * j + (i + 1) * i / 2 - 1));
      } 
    }
  }
  
  return psi_new;
}

/*** R
# Create artificial data
k <- 5
tt <- 4
n_covar <- (k - 1) * k / 2

# Constant parameters
psi <- matrix(1:(n_covar))
covar_vector_to_lower_triangular(psi, k, tt)

# Time varying parameters
psi <- matrix(1:(n_covar * tt))
covar_vector_to_lower_triangular(psi, k, tt)
*/