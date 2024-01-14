#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draws of Error Variances
//' 
//' Produces a draw of the constant diagonal error variance matrix of the
//' measurement equation of a state space model using an inverse gamma posterior density.
//' 
//' @param u a \eqn{KT \times 1} vector of errors.
//' @param shape_prior a \eqn{K \times 1} vector of prior shape parameters.
//' @param rate_prior a \eqn{K \times 1} vector of prior rate parameters.
//' @param inverse logical. If \code{TRUE}, the function returns the precision matrix,
//' i.e. the inverse of the variance matrix. Defaults to \code{FALSE}.
//' 
//' @details For a model with measurement equation
//' \deqn{y_t = Z_{t} a_t + u_t}
//' with \eqn{u_t \sim N(0, \Sigma_{u})}
//' the function produces a draw of the constant diagonal error variance matrix
//' \eqn{\Simga_u}.
//' 
//' @references
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
//' @examples
//' 
//' k <- 10 # Number of endogenous variables
//' tt <- 1000 # Number of observations
//' 
//' set.seed(1234) # Set RNG seed
//' 
//' # Generate artificial error series with N(0, 1)
//' u <- matrix(rnorm(k * tt))
//' 
//' # Define priors
//' shape_prior <- matrix(1, k)
//' rate_prior <- matrix(.0001, k)
//' 
//' # Obtain posterior draw
//' post_gamma_measurement_variance(u, shape_prior, rate_prior, inverse = FALSE)
//' 
//' @return A matrix.
//' 
// [[Rcpp::export]]
arma::sp_mat post_gamma_measurement_variance(const arma::vec u, const arma::vec shape_prior, const arma::vec rate_prior, const bool inverse) {
  
  // Get dimensions
  const int k = shape_prior.n_elem;
  const int tt = u.n_elem / shape_prior.n_elem;
  
  // Create vector of lagged values
  arma::vec a_v = arma::sum(arma::reshape(arma::square(u), k, tt), 1);
  const arma::vec scale_post = 1 / (rate_prior + a_v * 0.5);
  const arma::vec shape_post = shape_prior + tt * 0.5;
  arma::sp_mat result = arma::eye<arma::sp_mat>(k, k);
  for (int i = 0; i < k; i++) {
    // Inverse
    result(i, i) = arma::randg<double>(arma::distr_param(shape_post(i), scale_post(i)));
    if (!inverse) {
      result(i, i) = 1 / result(i, i);
    }
  }
  return result;
}

/*** R
k <- 10 # Number of artificial coefficients
tt <- 1000 # Number of observations

set.seed(1234) # Set RNG seed

# Generate artificial data according to a random walk
u <- matrix(rnorm(k * tt))

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
bla <- post_gamma_measurement_variance(u, shape_prior, rate_prior, inverse = FALSE)
bla

*/