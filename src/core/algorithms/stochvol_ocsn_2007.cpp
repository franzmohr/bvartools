// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "stochvol_ocsn_2007.h"

#include "bayests/arma.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

/**
 * @file stochvol_ocsn_2007.cpp
 * @brief Posterior draw of a stochastic volatility state with the normal
 *   mixture of Omori, Chib, Shephard and Nakajima (2007).
 */

namespace
{

// The ten-component normal mixture that approximates the log chi-squared
// distribution of log(u^2). It is what turns the non-linear measurement
// equation of a stochastic volatility model into a conditionally linear one.
const arma::rowvec::fixed<10> kMixtureWeight = {0.00609, 0.04775, 0.13057, 0.20674, 0.22715,
                                                0.18842, 0.12047, 0.05591, 0.01575, 0.00115};
const arma::rowvec::fixed<10> kMixtureMean = {1.92677, 1.34744, 0.73504, 0.02266, -0.85173,
                                              -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const arma::rowvec::fixed<10> kMixtureVariance = {0.11265, 0.17788, 0.26768, 0.40611, 0.62699,
                                                  0.98583, 1.57469, 2.54498, 4.16591, 7.33342};

const arma::uword kMixtureComponents = 10;

/// Everything the caller can get wrong about an argument is one message, so
/// that a failure names the argument rather than surfacing as an Armadillo
/// bounds error or a matrix of NaNs several draws later.
void require(bool ok, const std::string &what)
{
  if (!ok) {
    throw std::invalid_argument("stochvol_ocsn_2007: " + what);
  }
}

std::string dims(const arma::mat &m)
{
  return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

void require_length(const arma::vec &v, arma::uword n, const char *what)
{
  require(v.n_elem == n, std::string(what) + " must have " + std::to_string(n) +
                             " elements, got " + std::to_string(v.n_elem));
}

/// True if every element is a finite number, read off the exponent bits.
///
/// Not `arma::is_finite()` and not `std::isfinite()`: the library is compiled
/// with `-ffast-math`, which implies `-ffinite-math-only` and so licenses the
/// compiler to fold either of them to `true`. Armadillo says as much, with a
/// warning on stderr per call that a sampler would emit once per draw. Looking
/// at the bit pattern instead is an integer operation that no floating point
/// assumption can remove, and an all-ones exponent is exactly an infinity or a
/// NaN in IEEE-754 binary64.
bool all_finite(const arma::mat &m)
{
  static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

  const arma::uword n = m.n_elem;
  for (arma::uword i = 0; i < n; i++) {
    std::uint64_t bits;
    const double value = m[i];
    std::memcpy(&bits, &value, sizeof(bits));
    if ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL) {
      return false;
    }
  }
  return true;
}

/// Every element has to be finite and strictly positive: `sigma` is divided by
/// and `constant` is added under a logarithm, so a zero or a NaN here becomes
/// an infinity that only shows up as a broken draw further down.
void require_positive(const arma::vec &v, const char *what)
{
  require(all_finite(v), std::string(what) + " contains NaN or infinite values");
  require(v.min() > 0.0, std::string(what) + " must be strictly positive, but element " +
                             std::to_string(v.index_min() + 1) + " is " +
                             std::to_string(v.min()));
}

/// Draws the mixture indicator of every period for one variable.
///
/// The component probabilities are formed in logs and shifted by their row
/// maximum before they are exponentiated. The direct spelling -- weight times
/// `arma::normpdf`, then divide by the row sum -- underflows to ten zeros as
/// soon as an observation sits far enough out in the tails of all ten
/// components, and the normalisation then turns the row into NaNs and the draw
/// into an out-of-range index. After the shift the largest entry of every row
/// is exactly one, so the row sum can never be zero.
///
/// @param y_star T-vector of the transformed observations \f$\log(u_t^2 + c)\f$.
/// @param h_col T-vector of the current log-volatility.
/// @return T-vector of indices in `[0, 9]`.
arma::uvec draw_mixture_states(const arma::vec &y_star, const arma::vec &h_col)
{
  const arma::uword tt = y_star.n_elem;
  const arma::uword n = kMixtureComponents;

  // log p_j - 0.5 log sigma2_j - 0.5 (y*_t - h_t - mu_j)^2 / sigma2_j, up to
  // the constant that the normalisation removes anyway.
  arma::mat lq = arma::repmat(arma::log(kMixtureWeight) - 0.5 * arma::log(kMixtureVariance), tt, 1);
  lq -= 0.5 *
        arma::square(arma::repmat(y_star - h_col, 1, n) - arma::repmat(kMixtureMean, tt, 1)) /
        arma::repmat(kMixtureVariance, tt, 1);

  arma::mat q = lq.each_col() - arma::max(lq, 1);
  q = arma::exp(q);
  q = q.each_col() / arma::sum(q, 1);

  // Inverse c.d.f. draw: the number of components whose cumulated probability
  // still exceeds the uniform draw counts down from the component that is hit.
  arma::uvec s = n - arma::sum(arma::repmat(arma::randu<arma::vec>(tt), 1, n) <
                                   arma::cumsum(q, 1),
                               1);

  // A uniform draw above the last cumulated probability -- reachable by
  // rounding alone, since the row sums to one only to within floating point --
  // leaves no component hit and would index one past the last. Round it down to
  // the last component rather than letting `.elem()` throw.
  s.elem(arma::find(s >= n)).fill(n - 1);

  return s;
}

} // namespace

/**
 * @brief Draws the log-volatility of a stochastic volatility model.
 *
 * For a series of error terms \f$u_{it}\f$, \f$t = 1, \dots, T\f$, with
 * \f[
 *   u_{it} \sim N(0, \exp(h_{it})), \qquad
 *   h_{it} = h_{i,t-1} + v_{it}, \qquad v_{it} \sim N(0, \sigma_i),
 * \f]
 * the measurement equation is linearised by squaring and taking logarithms,
 * \f[
 *   \log(u_{it}^2 + c_i) = h_{it} + \varepsilon_{it},
 * \f]
 * where \f$\varepsilon_{it}\f$ is \f$\log \chi^2_1\f$ and \f$c_i\f$ is a small
 * offset that keeps the logarithm finite. That distribution is approximated by
 * the ten-component normal mixture of Omori, Chib, Shephard and Nakajima
 * (2007). Conditional on a mixture indicator per period the state space is
 * linear and Gaussian, so the whole path \f$h_{i1}, \dots, h_{iT}\f$ is drawn
 * in one block from its normal conditional posterior. Each column of `y` is
 * handled independently.
 *
 * @param y T x K matrix of error terms, one period per row. Must be finite;
 *   zeros are permitted because of the offset in `constant`.
 * @param h T x K matrix of the current log-volatility, used as the conditioning
 *   value of the mixture indicators. Must have the dimensions of `y`.
 * @param sigma K-vector of the current variances of the log-volatility
 *   innovations. Must be finite and strictly positive.
 * @param h_init K-vector of the initial states \f$h_{i0}\f$. Must be finite.
 * @param constant K-vector of the offsets \f$c_i\f$ added before the logarithm.
 *   Must be finite and strictly positive.
 *
 * @return T x K matrix of the drawn log-volatility. The argument `h` is left
 *   untouched, so a caller may pass the same matrix it assigns the result to.
 *
 * @throws std::invalid_argument if the arguments do not describe a stochastic
 *   volatility model of at least two periods and one variable, if their
 *   dimensions disagree, or if any of them contains a non-finite or
 *   non-positive value where the algorithm needs a finite or positive one. The
 *   message names the argument.
 * @throws std::runtime_error if a posterior precision turns out not to be
 *   positive definite, or if a drawn path is not finite. Either means the
 *   conditioning values have degenerated, and the message says which variable
 *   it was.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 *
 * @par Example
 * @code
 * arma::mat u = arma::randn<arma::mat>(100, 3);          // T = 100, K = 3
 * arma::mat h = arma::zeros<arma::mat>(100, 3);
 * arma::vec sigma = arma::vec(3).fill(0.05);
 * arma::vec h_init = arma::zeros<arma::vec>(3);
 * arma::vec constant = arma::vec(3).fill(1e-8);
 *
 * h = stochvol_ocsn_2007(u, h, sigma, h_init, constant);
 * @endcode
 */
arma::mat stochvol_ocsn_2007(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                             const arma::vec &h_init, const arma::vec &constant)
{
  // Checks
  require(y.n_rows == h.n_rows && y.n_cols == h.n_cols,
          "'h' must have the dimensions of 'y' (" + dims(y) + "), got " + dims(h));
  require(y.n_cols > 0, "'y' must have at least one column");
  require(y.n_rows > 1, "a stochastic volatility model needs at least two periods, got " +
                            std::to_string(y.n_rows));
  require(all_finite(y), "'y' contains NaN or infinite values");
  require(all_finite(h), "'h' contains NaN or infinite values");

  const arma::uword k = y.n_cols;
  const arma::uword tt = y.n_rows;

  require_length(sigma, k, "'sigma'");
  require_length(h_init, k, "'h_init'");
  require_length(constant, k, "'constant'");
  require_positive(sigma, "'sigma'");
  require_positive(constant, "'constant'");
  require(all_finite(h_init), "'h_init' contains NaN or infinite values");

  // The second difference operator of the random walk the log-volatility
  // follows, as a precision: hh = D'D with D the first difference. Sparse
  // because only three bands of the T x T matrix are ever populated.
  arma::sp_mat hh = arma::eye<arma::sp_mat>(tt, tt);
  hh.diag(-1) = -arma::ones<arma::vec>(tt - 1);
  hh = hh.t() * hh;

  arma::sp_mat sigs = arma::eye<arma::sp_mat>(tt, tt);
  const arma::vec vec_tt = arma::ones<arma::vec>(tt);

  arma::mat draw = h;

  for (arma::uword i = 0; i < k; i++) {
    const std::string variable = "variable " + std::to_string(i + 1);

    // Prepare series. Finite by construction: 'y' is finite and 'constant' is
    // strictly positive, so the argument of the logarithm is too.
    const arma::vec y_star = arma::log(arma::square(y.col(i)) + constant(i));

    // Sample the mixture indicators
    const arma::uvec s = draw_mixture_states(y_star, draw.col(i));

    // Sample the log-volatility in one block from N(V^-1 b, V^-1)
    const arma::sp_mat sigh_hh = hh / sigma(i);
    sigs.diag() = 1 / kMixtureVariance.elem(s);

    const arma::mat post_h_v(sigh_hh + sigs);
    const arma::vec rhs = sigh_hh * vec_tt * h_init(i) + sigs * (y_star - kMixtureMean.elem(s));

    // One Cholesky serves both the posterior mean and the draw: with
    // V = R'R the mean is R^-1 R^-T b and Cov(R^-1 z) = V^-1. The factorisation
    // is also the check -- V is positive definite by construction, being a
    // precision plus a positive diagonal, so a failure here means a
    // conditioning value has already degenerated.
    arma::mat r;
    if (!arma::chol(r, post_h_v)) {
      throw std::runtime_error("stochvol_ocsn_2007: the posterior precision of the log-volatility "
                               "of " + variable + " is not positive definite");
    }

    arma::vec forward, mean, noise;
    const bool solved = arma::solve(forward, arma::trimatl(r.t()), rhs) &&
                        arma::solve(mean, arma::trimatu(r), forward) &&
                        arma::solve(noise, arma::trimatu(r), arma::randn<arma::vec>(tt));
    if (!solved) {
      throw std::runtime_error("stochvol_ocsn_2007: the posterior of the log-volatility of " +
                               variable + " is singular");
    }

    const arma::vec h_i = mean + noise;
    if (!all_finite(h_i)) {
      throw std::runtime_error("stochvol_ocsn_2007: the drawn log-volatility of " + variable +
                               " is not finite");
    }

    draw.col(i) = h_i;
  }

  return draw;
}
