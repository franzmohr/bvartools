#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draws of Error Variances
//' 
//' Produces a draw of the constant diagonal error variance matrix of the
//' state equation of a state space model using an inverse gamma posterior density.
//' 
//' @param a a \eqn{KT \times 1} vector of time varying parameter draws.
//' @param a_init a \eqn{K \times 1} vector of initial states.
//' @param shape_prior a \eqn{K \times 1} vector of prior shape parameters.
//' @param rate_prior a \eqn{K \times 1} vector of prior rate parameters.
//' @param inverse logical. If \code{TRUE}, the function returns the precision matrix,
//' i.e. the inverse of the variance matrix. Defaults to \code{FALSE}.
//' 
//' @details For the state space model with state equation
//' \deqn{a_t = a_{t-1} + v}
//' and measurement equation
//' \deqn{y_t = Z_{t} a_t + u_t}
//' with \eqn{v_t \sim N(0, \Sigma_{v})} and \eqn{u_t \sim N(0, \Sigma_{u,t})}
//' the function produces a draw of the constant diagonal error variances matrix of the
//' state equation \eqn{\Simga_v}.
//' 
//' @references
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
//' @examples
//' k <- 10 # Number of artificial coefficients
//' tt <- 1000 # Number of observations
//' 
//' set.seed(1234) # Set RNG seed
//' 
//' # Generate artificial data according to a random walk
//' a <- matrix(rnorm(k), k, tt + 1)
//' for (i in 2:(tt + 1)) {
//'   a[, i] <- a[, i - 1] + rnorm(k, 0, sqrt(1 / 100))
//' }
//' 
//' a_init <- matrix(a[, 1]) # Define initial state
//' a <- matrix(a[, -1]) # Drop initial state from main sample and make vector
//' 
//' # Define priors
//' shape_prior <- matrix(1, k)
//' rate_prior <- matrix(.0001, k)
//' 
//' # Obtain posterior draw
//' post_gamma_state_variance(a, a_init, shape_prior, rate_prior, inverse = FALSE)
//' 
//' @return A matrix.
//' 
// [[Rcpp::export]]
arma::sp_mat post_gamma_state_variance(const arma::vec a, const arma::vec a_init, const arma::vec shape_prior, const arma::vec rate_prior, const bool inverse) {
  
  // Get dimensions
  const int k = a_init.n_elem;
  const int tt = a.n_elem / a_init.n_elem;
  
  // Create vector of lagged values
  arma::vec a_lag = a;
  a_lag.subvec(0, k - 1) = a_init;
  a_lag.subvec(k, k * tt - 1) = a.subvec(0, k * (tt - 1) - 1);
  arma::vec a_v = arma::sum(arma::reshape(arma::square(a - a_lag), k, tt), 1);
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
a <- matrix(rnorm(k), k, tt + 1)
for (i in 2:(tt + 1)) {
  a[, i] <- a[, i - 1] + rnorm(k, 0, sqrt(1 / 100))
}

a_init <- matrix(a[, 1]) # Define inital state
a <- a[, -1] # Drop initial state from main sample
a <- matrix(a)

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
bla <- post_gamma_state_variance(a, a_init, shape_prior, rate_prior, inverse = FALSE)
bla

*/