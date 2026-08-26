#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

#include "core/algorithms/kalman_durbin_koopman_2002.h"

//' Durbin and Koopman Simulation Smoother
//'
//' An implementation of the Kalman filter and backward smoothing
//' algorithm proposed by Durbin and Koopman (2002).
//'
//' @param y a \eqn{K \times T} matrix of endogenous variables.
//' @param z a \eqn{KT \times M} matrix of explanatory variables.
//' @param sigma_u the constant \eqn{K \times K} error variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{KT \times K} can be specified.
//' @param sigma_v the constant \eqn{M \times M} coefficient variance-covariance matrix.
//' For time varying variance-covariance matrices a \eqn{MT \times M} can be specified.
//' @param B the constant \eqn{M \times M} autocorrelation matrix of the transition
//' equation. For a time varying transition an \eqn{MT \times M} matrix can be specified.
//' @param a_init an M-dimensional vector of initial states.
//' @param P_init an \eqn{M \times M} variance-covariance matrix of the initial states.
//'
//' @details The function uses algorithm 2 from Durbin and Koopman (2002) to produce
//' a draw of the state vector \eqn{a_t} for \eqn{t = 1,...,T} for a state space model
//' with measurement equation
//' \deqn{y_t = Z_t a_t + u_t}
//' and transition equation
//' \deqn{a_{t + 1} = B_t a_{t} + v_t,}
//' where \eqn{u_t \sim N(0, \Sigma_{u,t})} and \eqn{v_t \sim N(0, \Sigma_{v,t})}.
//' \eqn{y_t} is a K-dimensional vector of endogenous variables and
//' \eqn{Z_t = z_t^{\prime} \otimes I_K} is a \eqn{K \times M} matrix of regressors with
//' \eqn{z_t} as a vector of regressors.
//'
//' The algorithm takes into account Jarociński (2015), where a possible missunderstanding
//' in the implementation of the algorithm of Durbin and Koopman (2002) is pointed out. Following
//' that note the function sets the mean of the initial state to zero in the first step of the algorithm.
//'
//' This is the routine the time varying parameter samplers of the package use. It
//' comes from the vendored BayesTS core, so there is one implementation of the
//' smoother rather than one for R and one for the samplers.
//'
//' The draw uses R's random number generator, so \code{\link{set.seed}} makes it
//' reproducible.
//'
//' @return A \eqn{M \times T+1} matrix of state vector draws. Column \eqn{i} is
//' the state the observation in column \eqn{i} of \code{y} loads on, for
//' \eqn{i = 1,...,T}, so a caller wanting one state per observation takes the
//' first \eqn{T} columns. The first column is not a state preceding the sample:
//' it is \eqn{a_1}, already conditioned on every observation. The last column is
//' the transition applied once past the end of the sample and is informed by no
//' observation; it is returned because the recursions build it on the way.
//'
//' @examples
//'
//' # Load data
//' data("e1")
//' data <- diff(log(e1))
//'
//' # Generate model data
//' temp <- create_bvarmodel(data = data, p = 2, deterministic = "const",
//'                          iterations = 1, burnin = 0)
//' y <- t(temp$data$train$y)
//' z <- temp$data$train$z
//' k <- nrow(y)
//' tt <- ncol(y)
//' m <- ncol(z)
//'
//' # Priors
//' a_mu_prior <- matrix(0, m)
//' a_v_i_prior <- diag(0.1, m)
//'
//' a_Q <- diag(.0001, m)
//'
//' # Initial value of Sigma
//' sigma <- tcrossprod(y) / tt
//' sigma_i <- solve(sigma)
//'
//' # Initial values for Kalman filter
//' y_init <- y * 0
//' a_filter <- matrix(0, m, tt + 1)
//'
//' # Initialise the Kalman filter
//' for (i in 1:tt) {
//'   y_init[, i] <- y[, i] - z[(i - 1) * k + 1:k,] %*% a_filter[, i]
//' }
//' a_init <- post_normal_sur(y = y_init, z = z, sigma_i = sigma_i,
//'                           a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
//' y_filter <- matrix(y) - z %*% a_init
//' y_filter <- matrix(y_filter, k) # Reshape
//'
//' # Kalman filter and backward smoother
//' a_filter <- kalman_durbin_koopman_2002(y = y_filter, z = z, sigma_u = sigma,
//'                                        sigma_v = a_Q, B = diag(1, m),
//'                                        a_init = matrix(0, m), P_init = a_Q)
//'
//' a <- a_filter + matrix(a_init, m, tt + 1)
//'
//' @references
//'
//' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for
//' state space time series analysis. \emph{Biometrika, 89}(3), 603--615.
//'
//' Jarociński, M. (2015). A note on implementing the Durbin and Koopman simulation
//' smoother. \emph{Computational Statistics and Data Analysis, 91}, 1--3.
//' \doi{10.1016/j.csda.2015.05.001}
//'
//' @export
// [[Rcpp::export(kalman_durbin_koopman_2002)]]
arma::mat kalman_durbin_koopman_2002_export(arma::mat y, arma::mat z,
                                            arma::mat sigma_u, arma::mat sigma_v,
                                            arma::mat B, arma::vec a_init, arma::mat P_init) {

  // The name this is exported under is the name of the function it calls, which
  // lives in the vendored core, so the two cannot both have it in C++. The
  // arguments are taken by value because the core takes four of them by
  // non-const reference, and because those copies are what the SEXP conversion
  // produces anyway.
  return kalman_durbin_koopman_2002(y, z, sigma_u, sigma_v, B, a_init, P_init);
}

/*** R

data("us_macrodata")

model <- create_bvarmodel(data = us_macrodata, p = 2,
                          deterministic = "const",
                          tvp = TRUE, iterations = 20, burnin = 10)

model <- add_priors(model,
                    coef = list(v_i = 1 / 9, v_i_det = 1 / 100, shape = 3, rate = .0001),
                    sigma = list(df = 3, scale = 1))

model <- add_initial_values(model)

k <- model[["model"]][["k"]]
y <- matrix(model[["data"]][["train"]][["y"]], k)
z <- model[["data"]][["train"]][["z"]]

sigma_u <- solve(model[["initial"]][["u_sigma_inv"]])
sigma_v <- as.matrix(solve(model[["initial"]][["a_v_inv"]]))
diag_b <- diag(1, ncol(z))
a_init <- model[["initial"]][["a_init"]]

kalman_durbin_koopman_2002(y, z, sigma_u, sigma_v, diag_b, a_init, sigma_v)

shape <- 3
rate <- .0001
set.seed(1); hist(rWishart(1000, df = shape, Sigma = 1 / matrix(rate)))
set.seed(1); hist(rgamma(1000, shape =  shape / 2, rate = matrix(rate) / 2))

*/
