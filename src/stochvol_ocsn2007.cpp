#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Stochastic Volatility
//' 
//' Produces a draw of log-volatilities based on Omori, Chib, Shephard and Nakajima (2007).
//' 
//' @param y a \eqn{T \times 1} vector containing the time series.
//' @param h a \eqn{T \times 1} vector of log-volatilities.
//' @param sigma a numeric of the variance of the log-volatilites.
//' @param h_init a numeric of the initial state of log-volatilities.
//' 
//' @details The function produces a posterior draw of the log-volatility \eqn{h} for the model
//' \deqn{y_{t} = e^{\frac{1}{2}h_t} \epsilon_{t},}
//' where \eqn{\epsilon_t \sim N(0, 1)} and \eqn{h_t} is assumed to evolve according to a random walk
//' \deqn{h_t = h_{t - 1} + u_t,}
//' with \eqn{u_t \sim N(0, \sigma^2)}.
//' 
//' The implementation follows the algorithm of Omori, Chib, Shephard and Nakajima (2007) and performs the
//' following steps:
//' \enumerate{
//'   \item Perform the transformation \eqn{y_t^* = ln(y_t^2 + 0.0001)}.
//'   \item Obtain a sample from the ten-component normal mixture for
//'   approximating the log-\eqn{\chi_1^2} distribution.
//'   \item Obtain a draw of log-volatilities.
//' }
//' 
//' @return A vector of log-volatility draws.
//' 
//' @examples
//' data("us_macrodata")
//' y <- us_macrodata[, "r"]
//' 
//' # Initialise log-volatilites
//' h_init <- log(var(y))
//' h <- rep(h_init, length(y))
//' 
//' # Obtain draw
//' stochvol_ocsn2007(y - mean(y), h, .05, h_init)
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
arma::vec stochvol_ocsn2007(arma::vec &y, arma::vec &h, double &sigma, double &h_init) {
  
  if (y.has_nan()) {
    Rcpp::stop("Argument 'y' contains NAs.");
  }
  
  // Prepare series
  y = log(arma::pow(y, 2) + 0.0001);
  arma::uword tt = y.n_elem;
  if (tt != h.n_elem) {
    Rcpp::stop("Arguments 'y' and 'h' do not have the same length.");
  }
  
  // Components of the mixture model
  arma::rowvec p_i(10), mu(10), sigma2(10);
  p_i(0) = 0.00609; mu(0) = 1.92677; sigma2(0) = 0.11265;
  p_i(1) = 0.04775; mu(1) = 1.34744; sigma2(1) = 0.17788;
  p_i(2) = 0.13057; mu(2) = 0.73504; sigma2(2) = 0.26768;
  p_i(3) = 0.20674; mu(3) = 0.02266; sigma2(3) = 0.40611;
  p_i(4) = 0.22715; mu(4) = -0.85173; sigma2(4) = 0.62699;
  p_i(5) = 0.18842; mu(5) = -1.97278; sigma2(5) = 0.98583;
  p_i(6) = 0.12047; mu(6) = -3.46788; sigma2(6) = 1.57469;
  p_i(7) = 0.05591; mu(7) = -5.55246; sigma2(7) = 2.54498;
  p_i(8) = 0.01575; mu(8) = -8.68384; sigma2(8) = 4.16591;
  p_i(9) = 0.00115; mu(9) = -14.65000; sigma2(9) = 7.33342;
  
  // Sample s
  arma::mat q = arma::repmat(p_i, tt, 1) % arma::normpdf(arma::repmat(y, 1, 10), arma::repmat(h, 1, 10) + arma::repmat(mu, tt, 1), arma::repmat(sqrt(sigma2), tt, 1));
  q = q / arma::repmat(sum(q, 1), 1, 10);
  arma::umat s = 10 - sum(arma::repmat(arma::randu<arma::vec>(tt), 1, 10) < cumsum(q, 1), 1);
  
  // Sample log-volatility
  arma::mat hh = arma::eye<arma::mat>(tt, tt);
  hh.diag(-1) = -arma::ones<arma::vec>(tt - 1);
  hh = hh.t() * hh;
  arma::mat sigh_hh = hh / sigma;
  arma::mat sigs = arma::eye<arma::mat>(tt, tt);
  sigs.diag() = 1 / sigma2.elem(s);
  arma::mat post_h_v = sigh_hh + sigs;
  arma::mat post_h_mu = arma::solve(post_h_v, sigh_hh * arma::ones<arma::vec>(tt) * h_init + sigs * (y - mu.elem(s)));
  
  return post_h_mu + arma::solve(arma::chol(arma::mat(post_h_v)), arma::randn<arma::vec>(tt));
}

/*** R

data("us_macrodata")
aud <- us_macrodata[, 1]
h_init <- log(var(aud))
h <- rep(h_init, length(aud))
stochvol_ocsn2007(aud - mean(aud), h, .05, h_init)

*/