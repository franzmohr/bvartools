#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draws of Error Variances
//' 
//' Produces a draw of the constant diagonal error variance matrix of the
//' state equation of a state space model using an inverse gamma posterior density.
//' 
//' @param a a \eqn{K \times T} matrix of time varying parameter draws.
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
//' a_init <- matrix(a[, 1]) # Define inital state
//' a <- a[, -1] # Drop initial state from main sample
//' 
//' # Define priors
//' shape_prior <- matrix(1, k)
//' rate_prior <- matrix(.0001, k)
//' 
//' # Obtain posterior draw
//' post_gamma_state_variance(a, a_init, shape_prior, rate_prior)
//' 
//' @return A matrix.
//' 
// [[Rcpp::export]]
arma::mat post_gamma_state_variance(arma::mat a, arma::vec a_init, arma::vec shape_prior, arma::vec rate_prior, bool inverse = false) {
  
  int k = a.n_rows;
  int tt = a.n_cols;
  arma::mat a_lag = a;
  a_lag.col(0) = a_init;
  a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
  arma::mat a_v = arma::trans(a - a_lag);
  arma::vec psi_sigma_v_post_scale = 1 / (rate_prior + arma::vectorise(arma::sum(arma::pow(a_v, 2))) * 0.5);
  arma::mat psi_sigma = arma::zeros<arma::mat>(k, k);
  arma::vec shape_post = shape_prior + tt * 0.5;
  for (int i = 0; i < k; i++) {
    psi_sigma(i, i) = arma::randg<double>(arma::distr_param(shape_post(i), psi_sigma_v_post_scale(i)));
    if (inverse) {
      psi_sigma(i, i) = 1 / psi_sigma(i, i);
    }
  }
  
  return psi_sigma;
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

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
post_gamma_state_variance(a, a_init, shape_prior, rate_prior, inverse = FALSE)

*/