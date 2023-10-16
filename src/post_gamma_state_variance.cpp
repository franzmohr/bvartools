#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Posterior Draws of Error Variances
//' 
//' Produces a draw of error variances from a gamma posterior density.
//' 
//' @param phi a \eqn{K \times T} matrix of time varying parameter draws.
//' @param phi_init a \eqn{K \times 1} vector of initial states.
//' @param shape_prior a \eqn{K \times 1} vector of prior shape parameters.
//' @param rate_prior a \eqn{K \times 1} vector of prior rate parameters.
//' 
//' @details The function produces a posterior draw of the variaces vector \eqn{a} for the model
//' Follow description in Chan eta al.
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
//' phi <- matrix(rnorm(k), k, tt + 1)
//' for (i in 2:(tt + 1)) {
//'   phi[, i] <- phi[, i - 1] + rnorm(k, 0, sqrt(1 / 100))
//' }
//' 
//' phi_init <- matrix(phi[, 1]) # Define inital state
//' phi <- phi[, -1] # Drop initial state from main sample
//' 
//' # Define priors
//' shape_prior <- matrix(1, k)
//' rate_prior <- matrix(.0001, k)
//' 
//' # Obtain posterior draw
//' post_gamma_state_variance(phi, phi_init, shape_prior, rate_prior)
//' 
//' @return A matrix.
//' 
// [[Rcpp::export]]
arma::mat post_gamma_state_variance(arma::mat phi, arma::vec phi_init, arma::vec shape_prior, arma::vec rate_prior) {
  
  int k = phi.n_rows;
  int tt = phi.n_cols;
  arma::mat phi_lag = phi;
  phi_lag.col(0) = phi_init;
  phi_lag.cols(1, tt - 1) = phi.cols(0, tt - 2);
  arma::mat phi_v = arma::trans(phi - phi_lag);
  arma::vec psi_sigma_v_post_scale = 1 / (rate_prior + arma::vectorise(arma::sum(arma::pow(phi_v, 2))) * 0.5);
  arma::mat psi_sigma_i = arma::zeros<arma::mat>(k, k);
  arma::vec shape_post = shape_prior + tt * 0.5;
  for (int i = 0; i < k; i++) {
    psi_sigma_i(i, i) = arma::randg<double>(arma::distr_param(shape_post(i), psi_sigma_v_post_scale(i)));
  }
  
  return psi_sigma_i;
}

/*** R
k <- 10 # Number of artificial coefficients
tt <- 1000 # Number of observations

set.seed(1234) # Set RNG seed

# Generate artificial data according to a random walk
phi <- matrix(rnorm(k), k, tt + 1)
for (i in 2:(tt + 1)) {
  phi[, i] <- phi[, i - 1] + rnorm(k, 0, sqrt(1 / 100))
}

phi_init <- matrix(phi[, 1]) # Define inital state
phi <- phi[, -1] # Drop initial state from main sample

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
post_gamma_state_variance(phi, phi_init, shape_prior, rate_prior)

*/