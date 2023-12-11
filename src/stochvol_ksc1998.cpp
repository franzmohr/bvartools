#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

//' Stochastic Volatility
//'
//' Produces a draw of log-volatilities.
//'
//' @param y a \eqn{T \times K} matrix containing the time series.
//' @param h a \eqn{T \times K} vector of log-volatilities.
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
//' The implementation is based on the algorithm of Kim, Shephard and Chip (1998) and performs the
//' following steps:
//' \enumerate{
//'   \item Perform the transformation \eqn{y_t^* = ln(y_t^2 + constant)}.
//'   \item Obtain a sample from the seven-component normal mixture for
//'   approximating the log-\eqn{\chi_1^2} distribution.
//'   \item Obtain a draw of log-volatilities.
//' }
//' 
//' The implementation follows the code provided on the website to the textbook
//' by Chan, Koop, Poirier, and Tobias (2019).
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
//' stochvol_ksc1998(y - mean(y), h, .05, h_init, 0.0001)
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
arma::mat stochvol_ksc1998(arma::mat y, arma::mat h, arma::vec sigma, arma::vec h_init, arma::vec constant) {
  
  if (y.has_nan()) {
    Rcpp::stop("Argument 'y' contains NAs.");
  }
  
  // Components of the mixture model
  arma::rowvec p_i(7), mu(7), sigma2(7);
  p_i(0) = 0.00730; p_i(1) = 0.10556; p_i(2) = 0.00002; p_i(3) = 0.04395; p_i(4) = 0.34001; p_i(5) = 0.24566; p_i(6) = 0.25750;
  mu(0) = -10.12999; mu(1) = -3.97281; mu(2) = -8.56686; mu(3) = 2.77786; mu(4) = 0.61942; mu(5) = 1.79518; mu(6) = -1.08819;
  mu = mu - 1.2704; // Already subtract the constant
  sigma2(0) = 5.79596; sigma2(1) = 2.61369; sigma2(2) = 5.17950; sigma2(3) = 0.16735; sigma2(4) = 0.64009; sigma2(5) = 0.34023; sigma2(6) = 1.26261;
  
  int k = y.n_cols;
  int tt = y.n_rows;
  arma::mat q, sigh_hh, sigs, post_h_v, post_h_mu;
  arma::umat s;
  
  arma::mat hh = arma::eye<arma::mat>(tt, tt);
  hh.diag(-1) = -arma::ones<arma::vec>(tt - 1);
  hh = hh.t() * hh;
  
  if (tt != h.n_rows) {
    Rcpp::stop("Arguments 'y' and 'h' do not have the same length.");
  }
  
  for (int i = 0; i < k; i++) {
    // Prepare series
    y.col(i) = log(arma::pow(y.col(i), 2) + constant(i));
    
    // Sample s
    q = arma::repmat(p_i, tt, 1) % arma::normpdf(arma::repmat(y.col(i), 1, 7), arma::repmat(h.col(i), 1, 7) + arma::repmat(mu, tt, 1), arma::repmat(sqrt(sigma2), tt, 1));
    q = q / arma::repmat(sum(q, 1), 1, 7);
    s = 7 - sum(arma::repmat(arma::randu<arma::vec>(tt), 1, 7) < cumsum(q, 1), 1);
    
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
stochvol_ksc1998(y, h, sigma_h, h_init, const)
***/
