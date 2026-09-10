// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "inverse_gaussian.h"
#include "bayests/arma.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

/**
 * @file inverse_gaussian.cpp
 * @brief Draws from the inverse Gaussian distribution.
 */

namespace
{

void require(bool ok, const std::string &what)
{
  if (!ok) {
    throw std::invalid_argument("inverse_gaussian: " + what);
  }
}

/// True if every element is a finite number, read off the exponent bits.
///
/// Not `arma::is_finite()` and not `std::isfinite()`, for the reason
/// `stochvol_mixture.h` sets out at length: a host that compiles these sources
/// with `-ffast-math` is licensed to fold either of them to `true`.
bool all_finite(const arma::vec &v)
{
  static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

  const arma::uword n = v.n_elem;
  for (arma::uword i = 0; i < n; i++) {
    std::uint64_t bits;
    const double value = v[i];
    std::memcpy(&bits, &value, sizeof(bits));
    if ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL) {
      return false;
    }
  }
  return true;
}

} // namespace

/**
 * @brief Draws from the inverse Gaussian distribution, one draw per element.
 *
 * For each \f$i\f$ draws \f$x_i \sim IG(\mu_i, \lambda_i)\f$ with density
 * \f[
 *   f(x) = \sqrt{\frac{\lambda}{2 \pi x^3}}
 *          \exp\left(-\frac{\lambda (x - \mu)^2}{2 \mu^2 x}\right),
 *          \qquad x > 0,
 * \f]
 * so \f$E[x] = \mu\f$ and \f$Var[x] = \mu^3 / \lambda\f$.
 *
 * The transformation of Michael, Schucany and Haas (1976): a chi-square variate
 * gives two roots of a quadratic, one of which is selected by a uniform. Exact
 * rather than a rejection scheme, so the number of random variates consumed per
 * draw is fixed -- which is what keeps a chain reproducible under a fixed seed
 * regardless of the values passed in.
 *
 * Written for the latent scale of an asymmetric Laplace model. There the full
 * conditional of the scale \f$w\f$ is generalised inverse Gaussian with index
 * \f$1/2\f$, which is the *reciprocal* of an inverse Gaussian -- so a caller
 * draws from here and inverts. It is the reciprocal rather than the variate
 * itself because the index is positive; at \f$-1/2\f$ it would be the variate.
 *
 * @param mu the means, all strictly positive.
 * @param lambda the shape parameters, all strictly positive. Either the same
 *   length as `mu`, or a single element that applies to every draw.
 *
 * @return a vector of draws the length of `mu`, every element strictly positive.
 *
 * @throws std::invalid_argument if either argument is empty, holds a
 *   non-positive or non-finite value, or if `lambda` is neither a scalar nor the
 *   length of `mu`.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 *
 * Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random
 * variates using transformations with multiple roots. The American
 * Statistician, 30(2), 88-90.
 */
arma::vec inverse_gaussian(const arma::vec &mu, const arma::vec &lambda)
{
  const arma::uword n = mu.n_elem;
  require(n > 0, "'mu' is empty");
  require(lambda.n_elem == n || lambda.n_elem == 1,
          "'lambda' must have 1 or " + std::to_string(n) + " elements, got " +
              std::to_string(lambda.n_elem));
  require(all_finite(mu), "'mu' contains NaN or infinite values");
  require(all_finite(lambda), "'lambda' contains NaN or infinite values");
  require(mu.min() > 0.0, "'mu' must be strictly positive");
  require(lambda.min() > 0.0, "'lambda' must be strictly positive");

  const bool lambda_is_scalar = lambda.n_elem == 1;

  arma::vec out(n);
  for (arma::uword i = 0; i < n; i++) {
    const double m = mu[i];
    const double l = lambda_is_scalar ? lambda[0] : lambda[i];

    // A chi-square variate on one degree of freedom, which is the square of a
    // standard normal.
    const double normal = arma::randn<double>();
    const double v = normal * normal;

    // The smaller root of the quadratic the transformation inverts. The form
    // usually written down,
    //
    //     m + m^2 v / (2l) - (m / (2l)) sqrt(4 m l v + m^2 v^2),
    //
    // subtracts two nearly equal quantities once m*v is large and can land at
    // or below zero, which is outside the support. Multiplying through by the
    // conjugate gives the algebraically identical
    //
    //     4 m^2 l v / (m v + sqrt(m^2 v^2 + 4 m l v))^2,
    //
    // in which every term is positive, so no root can be lost to cancellation
    // and no fallback is needed.
    double x;
    if (v > 0.0) {
      const double s = std::sqrt(m * m * v * v + 4.0 * m * l * v);
      const double denominator = m * v + s;
      x = 4.0 * m * m * l * v / (denominator * denominator);
    } else {
      // The limit of the expression above, and the value the transformation
      // takes when the normal draw is exactly zero.
      x = m;
    }

    const double u = arma::randu<double>();
    out[i] = (u <= m / (m + x)) ? x : m * m / x;
  }

  return out;
}
