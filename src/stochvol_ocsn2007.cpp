#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Stochastic Volatility
//' 
//' Produces a draw of log-volatilities based on Omori, Chib, Shephard and Nakajima (2007).
//' 
//' @param y a \eqn{T \times K} matrix containing the time series.
//' @param h a \eqn{T \times K} vector of the current draw of log-volatilities.
//' @param sigma a \eqn{K \times 1} vector of variances of log-volatilities,
//' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
//' @param h_init a \eqn{K \times 1} vector of the initial states of log-volatilities,
//' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
//' @param constant a \eqn{K \times 1} vector of constants that should be added to \eqn{y^2}
//' before taking the natural logarithm. The \eqn{i}th element corresponds to
//' the \eqn{i}th column in \code{y}. See 'Details'.
//' 
//' @details For each column in \code{y} the function produces a posterior
//' draw of the log-volatility \eqn{h} for the model
//' \deqn{y_{t} = e^{\frac{1}{2}h_t} \epsilon_{t},}
//' where \eqn{\epsilon_t \sim N(0, 1)} and \eqn{h_t} is assumed to evolve according to a random walk
//' \deqn{h_t = h_{t - 1} + u_t,}
//' with \eqn{u_t \sim N(0, \sigma^2)}.
//' 
//' The implementation follows the algorithm of Omori, Chib, Shephard and Nakajima (2007) and performs the
//' following steps:
//' \enumerate{
//'   \item Perform the transformation \eqn{y_t^* = ln(y_t^2 + constant)}.
//'   \item Obtain a sample from the ten-component normal mixture for
//'   approximating the log-\eqn{\chi_1^2} distribution.
//'   \item Obtain a draw of log-volatilities.
//' }
//' 
//' The implementation is an adaption of the code provided on the website to the textbook
//' by Chan, Koop, Poirier, and Tobias (2019).
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
//' stochvol_ocsn2007(y - mean(y), h, matrix(.05), h_init, matrix(0.0001))
//' 
//' @references
//' 
//' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
//' (2nd ed.). Cambridge: Cambridge University Press.
//' 
//' Omori, Y., Chib, S., Shephard, N., & Nakajima, J. (2007). Stochastic volatiltiy with leverage. Fast and efficient likelihood inference.
//' \emph{Journal of Econometrics 140}(2), 425--449. \doi{10.1016/j.jeconom.2006.07.008}
//' 
// [[Rcpp::export]]
arma::mat stochvol_ocsn2007(arma::mat y, arma::mat h, arma::vec sigma, arma::vec h_init, arma::vec constant) {
  
  // Checks
  if (y.has_nan()) {
    Rcpp::stop("Argument 'y' contains NAs.");
  }
  if (y.n_rows != h.n_rows) {
    Rcpp::stop("Arguments 'y' and 'h' do not have the same number of rows.");
  }
  if (y.n_cols != h.n_cols) {
    Rcpp::stop("Arguments 'y' and 'h' do not have the same number of columns.");
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
  
  int k = y.n_cols;
  int tt = y.n_rows;
  arma::mat q, sigh_hh, sigs, post_h_v, post_h_mu;
  arma::umat s;
  
  arma::mat hh = arma::eye<arma::mat>(tt, tt);
  hh.diag(-1) = -arma::ones<arma::vec>(tt - 1);
  hh = hh.t() * hh;
  
  for (int i = 0; i < k; i++) {
    // Prepare series
    y.col(i) = log(arma::pow(y.col(i), 2) + constant(i));
    
    // Sample s
    q = arma::repmat(p_i, tt, 1) % arma::normpdf(arma::repmat(y.col(i), 1, 10), arma::repmat(h.col(i), 1, 10) + arma::repmat(mu, tt, 1), arma::repmat(sqrt(sigma2), tt, 1));
    q = q / arma::repmat(sum(q, 1), 1, 10);
    s = 10 - sum(arma::repmat(arma::randu<arma::vec>(tt), 1, 10) < cumsum(q, 1), 1);
    
    // Sample log-volatility
    sigh_hh = hh / sigma(i);
    sigs = arma::eye<arma::mat>(tt, tt);
    sigs.diag() = 1 / sigma2.elem(s);
    post_h_v = sigh_hh + sigs;
    post_h_mu = arma::solve(post_h_v, sigh_hh * arma::ones<arma::vec>(tt) * h_init(i) + sigs * (y.col(i) - mu.elem(s)));
    h.col(i) = post_h_mu + arma::solve(arma::chol(post_h_v), arma::randn<arma::vec>(tt));
  }
  
  return h;
}

/*** R

data("us_macrodata")
y <- us_macrodata
h_init <- log(diag(var(y)))
h <- t(matrix(h_init, 3, nrow(y)))
sigma_h <- rep(.05, 3)
const <- rep(.0001, 3)
stochvol_ocsn2007(y, h, sigma_h, h_init, const)

***/
