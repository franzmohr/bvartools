// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// Pareto smoothed importance sampling, Vehtari, Simpson, Gelman, Yao and
// Gabry (2024, Journal of Machine Learning Research 25(72), 1-58).
//
// WHAT IT IS FOR HERE. add_sign_zero_restrictions() reweights each draw by the
// ratio of two volume elements, and that ratio is unbounded: nothing in the
// algorithm stops one draw from landing where the proposal put almost no
// probability and the target a great deal. When that happens the weight of that
// draw is enormous, the effective sample size collapses, and the resample is a
// few copies of one draw wearing the clothes of a posterior. It is not a rare
// accident -- a five-variable Austrian VAR under five zero restrictions produced
// an effective sample size of 4 out of 487 accepted draws, one of them holding
// 46 per cent of the total weight.
//
// The importance weights of such a sampler typically follow a generalised
// Pareto distribution in their upper tail. PSIS fits one to the largest few
// hundred weights and replaces them by the quantiles of the fit. The fitted
// shape parameter k is what makes this more than a smoothing trick: it is an
// estimate of how heavy that tail is, and it says when the sample cannot be
// trusted at all rather than leaving the user to guess from a small effective
// sample size. A finite variance requires k < 1/2, and the practical threshold
// the authors recommend is 0.7, above which neither the estimate nor its
// effective sample size means very much.
//
// WHAT IT COSTS. The smoothed estimator is biased where the raw one is not.
// The bias is small and the variance reduction is large, which is the trade the
// paper argues for and the one taken here -- but it is a trade, which is why
// add_sign_zero_restrictions() can be asked for the raw weights of the original
// Algorithm 3 instead.
//
// This is a translation of the reference implementation in the loo package,
// which the authors maintain, rather than an independent derivation. It lives
// here rather than as a dependency because it is short, because reweighting an
// already estimated model is this package's own business, and because a student
// reading why their draws were reweighted should find the answer in the source
// rather than in another library.

namespace {

double log_sum_exp(const std::vector<double> &x) {
  double max_value = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < x.size(); i++) {
    if (x[i] > max_value) {
      max_value = x[i];
    }
  }
  if (!std::isfinite(max_value)) {
    return max_value;
  }
  double total = 0.0;
  for (std::size_t i = 0; i < x.size(); i++) {
    total += std::exp(x[i] - max_value);
  }
  return max_value + std::log(total);
}

// The profile log-likelihood of the generalised Pareto distribution at one
// value of theta, divided by the sample size. `x` is sorted and positive.
double profile_loglik(double theta, const std::vector<double> &x) {
  const double a = -theta;
  double k = 0.0;
  for (std::size_t i = 0; i < x.size(); i++) {
    k += std::log1p(a * x[i]);
  }
  k /= static_cast<double>(x.size());
  // theta runs over both signs on this grid, and a and k always carry the same
  // one, so what has to be guarded is not the sign of either but the ratio
  // being a number the logarithm accepts.
  const double ratio = a / k;
  if (!std::isfinite(ratio) || !(ratio > 0.0)) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::log(ratio) - k - 1.0;
}

struct GpdFit {
  double k;
  double sigma;
  bool ok;
};

// Zhang and Stephens (2009), the empirical-Bayes estimator the reference
// implementation uses: a grid of candidate theta values, each weighted by its
// own profile likelihood, averaged rather than maximised. `x` must be sorted
// ascending and strictly positive.
GpdFit fit_generalised_pareto(const std::vector<double> &x) {
  GpdFit out;
  out.k = std::numeric_limits<double>::quiet_NaN();
  out.sigma = std::numeric_limits<double>::quiet_NaN();
  out.ok = false;

  const std::size_t n = x.size();
  if (n < 5 || !(x.back() > 0.0)) {
    return out;
  }

  // The quarter point of the sample, which sets the scale of the grid.
  const std::size_t quarter =
    static_cast<std::size_t>(std::floor(static_cast<double>(n) / 4.0 + 0.5));
  const double x_star = x[quarter > 0 ? quarter - 1 : 0];
  if (!(x_star > 0.0)) {
    return out;
  }

  const int grid = 30 + static_cast<int>(std::floor(std::sqrt(static_cast<double>(n))));
  std::vector<double> theta(grid);
  std::vector<double> log_lik(grid);
  for (int j = 0; j < grid; j++) {
    const double jj = static_cast<double>(j) + 1.0;
    theta[j] = 1.0 / x.back() +
      (1.0 - std::sqrt(static_cast<double>(grid) / (jj - 0.5))) / (3.0 * x_star);
    log_lik[j] = static_cast<double>(n) * profile_loglik(theta[j], x);
  }

  const double normaliser = log_sum_exp(log_lik);
  if (!std::isfinite(normaliser)) {
    return out;
  }

  double theta_hat = 0.0;
  for (int j = 0; j < grid; j++) {
    theta_hat += theta[j] * std::exp(log_lik[j] - normaliser);
  }

  double k = 0.0;
  for (std::size_t i = 0; i < n; i++) {
    k += std::log1p(-theta_hat * x[i]);
  }
  k /= static_cast<double>(n);
  const double sigma = -k / theta_hat;

  // The weakly informative prior on k of the reference implementation, which
  // pulls a shape estimated from few tail draws towards 0.5 rather than letting
  // it run away on the strength of a handful of points.
  const double a = 10.0;
  const double n_plus_a = static_cast<double>(n) + a;
  out.k = k * static_cast<double>(n) / n_plus_a + a * 0.5 / n_plus_a;
  out.sigma = sigma;
  out.ok = std::isfinite(out.k) && std::isfinite(sigma) && sigma > 0.0;
  return out;
}

double generalised_pareto_quantile(double p, double k, double sigma) {
  if (std::abs(k) < 1e-30) {
    return -sigma * std::log1p(-p);
  }
  return sigma * std::expm1(-k * std::log1p(-p)) / k;
}

}  // namespace

//' Pareto smoothed importance weights
//'
//' @param log_weight numeric. The unnormalised log importance weights. Draws
//'   that carry no weight are passed as \code{-Inf} and are neither smoothed nor
//'   counted towards the tail.
//'
//' @return A list with \code{log_weights}, the smoothed and normalised log
//'   weights in the order they were given, and \code{pareto_k}, the shape of the
//'   generalised Pareto distribution fitted to the tail. \code{pareto_k} is
//'   \code{NA} when there were too few finite weights to fit one, in which case
//'   the weights come back normalised but unsmoothed.
//'
//' @noRd
// [[Rcpp::export(.psis_smooth)]]
Rcpp::List psis_smooth(Rcpp::NumericVector log_weight) {

  const int total = log_weight.size();
  Rcpp::NumericVector out(total, R_NegInf);
  double pareto_k = NA_REAL;

  // Only the draws that carry weight take part. A rejected draw is not a draw
  // with a very small weight -- it is not in the sample at all, and letting it
  // into the tail fit would describe a distribution that has no draws in it.
  std::vector<int> keep;
  keep.reserve(total);
  for (int i = 0; i < total; i++) {
    if (R_finite(log_weight[i])) {
      keep.push_back(i);
    }
  }

  const int n = static_cast<int>(keep.size());
  if (n == 0) {
    return Rcpp::List::create(Rcpp::Named("log_weights") = out,
                              Rcpp::Named("pareto_k") = pareto_k);
  }

  std::vector<double> lw(n);
  for (int i = 0; i < n; i++) {
    lw[i] = log_weight[keep[i]];
  }

  // Everything below is on the scale of the largest weight, which keeps exp()
  // of a log weight in the hundreds from overflowing.
  const double max_lw = *std::max_element(lw.begin(), lw.end());
  for (int i = 0; i < n; i++) {
    lw[i] -= max_lw;
  }

  // The tail the reference implementation smooths: a fifth of the sample, or
  // three times its square root, whichever is smaller.
  const int tail_length = static_cast<int>(std::ceil(
    std::min(0.2 * static_cast<double>(n), 3.0 * std::sqrt(static_cast<double>(n)))));

  // Fewer than five tail draws cannot say anything about the shape of a tail.
  if (tail_length >= 5 && n > tail_length) {
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&lw](int a, int b) { return lw[a] < lw[b]; });

    const int first_tail = n - tail_length;
    const double cutoff = lw[order[first_tail - 1]];

    // Exceedances over the cutoff, ascending by construction of `order`.
    const double exp_cutoff = std::exp(cutoff);
    std::vector<double> excess(tail_length);
    for (int j = 0; j < tail_length; j++) {
      excess[j] = std::exp(lw[order[first_tail + j]]) - exp_cutoff;
    }

    // A tail whose draws are all the same value has no shape to fit, and the
    // fit would divide by zero rather than fail loudly.
    const bool degenerate = (excess.back() - excess.front()) <
      std::numeric_limits<double>::epsilon() * 100.0;

    if (!degenerate) {
      const GpdFit fit = fit_generalised_pareto(excess);
      if (fit.ok) {
        pareto_k = fit.k;
        for (int j = 0; j < tail_length; j++) {
          const double p = (static_cast<double>(j) + 0.5) /
            static_cast<double>(tail_length);
          const double q =
            generalised_pareto_quantile(p, fit.k, fit.sigma) + exp_cutoff;
          if (q > 0.0) {
            lw[order[first_tail + j]] = std::log(q);
          }
        }
      }
    }
  }

  // No smoothed weight may exceed the largest weight actually drawn: the fitted
  // quantiles extrapolate beyond the sample and the truncation is what keeps
  // that extrapolation from inventing a draw heavier than any that occurred.
  for (int i = 0; i < n; i++) {
    if (lw[i] > 0.0) {
      lw[i] = 0.0;
    }
  }

  const double normaliser = log_sum_exp(lw);
  for (int i = 0; i < n; i++) {
    out[keep[i]] = lw[i] - normaliser;
  }

  return Rcpp::List::create(Rcpp::Named("log_weights") = out,
                            Rcpp::Named("pareto_k") = pareto_k);
}
