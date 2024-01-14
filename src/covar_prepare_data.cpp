#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Covariance: Data Preparation
//'
//' Convenience function, which generates the input data for posterior
//' simulation of covariance parameters.
//'
//' @param y a \eqn{KT \times 1} vector of input data.
//' @param omega_i a \eqn{K \times K} or \eqn{KT \times KT} matrix of error variances.
//' The matrix must be sparse.
//' @param k an integer of the number of endogenous variables.
//' @param tt an integer of the number of observations.
//' @param tvp logical indicating if the SUR matrix with the values of regressors
//' should be prepared for the estimation of constant or time varying parameters.
//' 
//' @details For the model
//' \deqn{y_t = Z_{t} a_t + u_t}
//' with \eqn{u_t \sim N(0, \Psi \Omega_{t} \Psi^{\prime})} and \eqn{\Omega_{t}}
//' as a diagonal matrix of error variances, the function produces
//' the input data for the posterior simulation of the lower triangular covariance coefficients
//' of \eqn{\Psi} as presented in Primiceri (2005).
//' 
//' @return A list with three elements:
//' \item{y}{The prepared vector of endogenous variables.}
//' \item{z}{The prepared matrix of regressors.}
//' \item{omega_i}{The prepared diagonal matrix of measurement error variances.}
//' All matrices are returned as sparse matrices.
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
//' k <- 3
//' tt <- 4
//' u <- matrix(1:(k * tt))
//' omega_i <- Matrix(diag(1:3, k))
//' 
//' # Generate input data (constant parameters)
//' covar_prepare_data(u, omega_i, k, tt, FALSE)
//' 
//' # Generate input data (time varying parameters)
//' covar_prepare_data(u, omega_i, k, tt, TRUE)
//' 
// [[Rcpp::export(covar_prepare_data)]]
Rcpp::List covar_prepare_data(const arma::vec y, const arma::sp_mat omega_i, const arma::uword k, const int tt, const bool tvp) {
  
  if (k < 2) {
    Rcpp::stop("Argument 'k' must be larger than 1.");
  }
  
  const int n_covar = k * (k - 1) / 2;
  
  arma::sp_mat omega_i_new, z;
  
  if (tvp) {
    z = arma::zeros<arma::sp_mat>((k - 1) * tt, n_covar * tt);
    for (int j = 0; j < tt; j++) {
      for (arma::uword i = 1; i < k; i++) {
        z.submat(j * (k - 1) + i - 1,
                 j * n_covar + i * (i - 1) / 2,
                 j * (k - 1) + i - 1,
                 j * n_covar + (i + 1) * i / 2 - 1) = -arma::trans(y.subvec(j * k, j * k + i - 1));
      }
    }
  } else {
    z = arma::zeros<arma::sp_mat>((k - 1) * tt, n_covar);
    for (int j = 0; j < tt; j++) {
      for (arma::uword i = 1; i < k; i++) {
        z.submat(j * (k - 1) + i - 1,
                 i * (i - 1) / 2,
                 j * (k - 1) + i - 1,
                 (i + 1) * i / 2 - 1) = -arma::trans(y.subvec(j * k, j * k + i - 1));
      }
    }
  }
  
  if (omega_i.n_cols == k) {
    omega_i_new = arma::kron(arma::eye<arma::sp_mat>(tt, tt), omega_i.submat(1, 1, k - 1, k - 1));
  } else {
    if (omega_i.n_cols == k * tt) {
      omega_i_new = arma::zeros<arma::sp_mat>((k - 1) * tt, (k - 1) * tt);
      for (int j = 0; j < tt; j++) {
        omega_i_new.submat(j * (k - 1),
                           j * (k - 1),
                           (j + 1) * (k - 1) - 1,
                           (j + 1) * (k - 1) - 1) = omega_i.submat(j * k + 1, j * k + 1, (j + 1) * k - 1, (j + 1) * k - 1);
      }
    } 
  }
  
  arma::mat y_new = arma::reshape(y, k, tt);
  
  return Rcpp::List::create(Rcpp::Named("y") = arma::vectorise(y_new.rows(1, k - 1)),
                            Rcpp::Named("z") = z,
                            Rcpp::Named("omega_i") = omega_i_new);
}

/*** R
k <- 3
tt <- 4
u <- matrix(1:(k * tt))

omega_i <- Matrix(diag(1:k, k))

prepare_covar_data(u, omega_i, k, tt, FALSE)

prepare_covar_data(u, omega_i, k, tt, TRUE)

# Time varying errors
new_omega_i <- Matrix(0, k * tt, k * tt)
for (i in 1:tt) {
  new_omega_i[(i - 1) * k + 1:k, (i - 1) * k + 1:k] <- omega_i
}

prepare_covar_data(u, new_omega_i, k, tt, FALSE)

prepare_covar_data(u, new_omega_i, k, tt, TRUE)
*/