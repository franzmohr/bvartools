#include "../inst/include/bvartools.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Stochastic Volatility
//' 
//' Produces a draw of log-volatilities.
//' 
//' @param y a \eqn{T \times 1} vector containing the time series.
//' @param h a \eqn{T \times 1} vector of log-volatilities.
//' @param sigma a numeric of the variance of the log-volatilites.
//' @param h_init a numeric of the initial state of log-volatilities.
//' @param constant a numeric of the constant that should be added to \eqn{y^2}
//' before taking the natural logarithm.
//' 
//' @details The function is a wrapper for function \code{\link{stochvol_ksc1998}}.
//' 
//' @return A vector of log-volatility draws.
//' 
//' @examples
//' data("us_macrodata")
//' y <- matrix(us_macrodata[, "r"])
//' 
//' # Initialise log-volatilites
//' h_init <- matrix(log(var(y)))
//' h <- matrix(rep(h_init, length(y)))
//' 
//' # Obtain draw
//' stoch_vol(y - mean(y), h, matrix(.05), h_init, matrix(0.0001))
//' 
//' @references
//' 
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
//' Kim, S., Shephard, N., & Chib, S. (1998). Stochastic volatility. Likelihood inference and comparison
//' with ARCH models. \emph{Review of Economic Studies 65}(3), 361--393. \doi{10.1111/1467-937X.00050}
//' 
// [[Rcpp::export]]
arma::vec stoch_vol(arma::mat y, arma::mat h, arma::vec sigma, arma::vec h_init, arma::vec constant) {
  return bvartools::stochvol_ksc1998(y, h, sigma, h_init, constant);
}

/*** R

data("us_macrodata")
y <- us_macrodata[, 1]
h_init <- log(var(y))
h <- rep(h_init, length(y))
stoch_vol(y - mean(y), h, .05, h_init, 0.0001)

***/