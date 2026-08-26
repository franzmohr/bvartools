#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

#include "core/algorithms/stochvol_ocsn_2007.h"

//' Stochastic Volatility
//'
//' Produces a draw of log-volatilities based on Omori, Chib, Shephard and Nakajima (2007).
//'
//' @param y a \eqn{T \times K} matrix containing the time series.
//' @param h a \eqn{T \times K} matrix of the current draw of log-volatilities.
//' @param sigma a \eqn{K \times 1} vector of variances of log-volatilities,
//' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
//' Must be finite and positive.
//' @param h_init a \eqn{K \times 1} vector of the initial states of log-volatilities,
//' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
//' @param constant a \eqn{K \times 1} vector of constants that should be added to \eqn{y^2}
//' before taking the natural logarithm. The \eqn{i}th element corresponds to
//' the \eqn{i}th column in \code{y}. Must be finite and positive. See 'Details'.
//'
//' @details For each column in \code{y} the function produces a posterior
//' draw of the log-volatility \eqn{h} for the model
//' \deqn{y_{t} = e^{\frac{1}{2}h_t} \epsilon_{t},}
//' where \eqn{\epsilon_t \sim N(0, 1)} and \eqn{h_t} is assumed to evolve according to a random walk
//' \deqn{h_t = h_{t - 1} + u_t,}
//' with \eqn{u_t \sim N(0, \sigma^2)}.
//'
//' The implementation follows the algorithm of Omori, Chib, Shephard and Nakajima (2007)
//' and performs the following steps:
//' \enumerate{
//'   \item Perform the transformation \eqn{y_t^* = ln(y_t^2 + constant)}.
//'   \item Obtain a sample from the ten-component normal mixture for
//'   approximating the log-\eqn{\chi_1^2} distribution.
//'   \item Obtain a draw of log-volatilities.
//' }
//'
//' The ten components approximate the log-\eqn{\chi_1^2} distribution more
//' closely than the seven of \code{\link{stochvol_ksc_1998}}, at a proportional
//' cost in the first step.
//'
//' The posterior precision of the log-volatility is tridiagonal -- the random
//' walk contributes the two bands of \eqn{D^\prime D} and the mixture only a
//' diagonal -- so the draw is obtained from a banded Cholesky factorisation,
//' which costs \eqn{O(T)} rather than the \eqn{O(T^3)} of a dense one.
//'
//' The probabilities of the mixture components in step 2 are obtained in logs and
//' shifted by their period-wise maximum before they are exponentiated. Weighting
//' the component densities directly underflows to zero for all ten components
//' once an observation lies far enough out in their tails, which leaves the
//' probabilities of that period undefined and selects a component that does not
//' exist.
//'
//' The arguments are checked before the first draw and a violation raises an
//' error. \code{sigma} and \code{constant} have to be finite and positive, since
//' the algorithm divides by the former and takes the logarithm of a sum
//' containing the latter, and all three of \code{sigma}, \code{h_init} and
//' \code{constant} have to have one element per column of \code{y}.
//'
//' The draw uses R's random number generator, so \code{\link{set.seed}} makes it
//' reproducible.
//'
//' @return A \eqn{T \times K} matrix of log-volatility draws.
//'
//' @references
//'
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//'
//' Kim, S., Shephard, N., & Chib, S. (1998). Stochastic volatility. Likelihood inference and comparison
//' with ARCH models. \emph{Review of Economic Studies 65}(3), 361--393. \doi{10.1111/1467-937X.00050}
//'
//' Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic volatiltiy with leverage. Fast and efficient likelihood inference.
//' \emph{Journal of Econometrics 140}(2), 425--449. \doi{10.1016/j.jeconom.2006.07.008}
//'
//' @examples
//' data("us_macrodata")
//' y <- diff(us_macrodata)
//' h_init <- log(diag(var(y)))
//' h <- t(matrix(h_init, 3, nrow(y)))
//' sigma_h <- rep(.05, 3)
//' const <- rep(.0001, 3)
//' stochvol_ocsn_2007(y, h, sigma_h, h_init, const)
//'
//' @export
// [[Rcpp::export(stochvol_ocsn_2007)]]
arma::mat stochvol_ocsn_2007_export(const arma::mat &y, const arma::mat &h,
                                    const arma::vec &sigma, const arma::vec &h_init,
                                    const arma::vec &constant) {

  // The name this is exported under is the name of the function it calls, which
  // lives in the vendored core, so the two cannot both have it in C++.
  return stochvol_ocsn_2007(y, h, sigma, h_init, constant);
}
