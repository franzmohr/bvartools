// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_ALD_SUPPORT_H
#define BAYESTS_CORE_MODELS_ALD_SUPPORT_H

#include "bayests/arma.h"
#include "core/algorithms/inverse_gaussian.h"

#include <cmath>
#include <stdexcept>

/// @file ald_support.h
/// @brief The blocks a quantile regression model draws that a mean regression
///        does not, shared by VarNormalAld and VarTvpAld.
///
/// Everything here follows from one representation. For quantile q, with
///
///     theta = (1 - 2q) / (q (1 - q)),   tau2 = 2 / (q (1 - q)),
///
/// the asymmetric Laplace with scale s is the marginal of
///
///     u | w ~ N(theta w, tau2 s w),   w ~ Exp(1 / s),
///
/// so a model whose errors are asymmetric Laplace is, conditional on w, a
/// weighted normal model -- and the two draws below are all that separates the
/// two. Kozumi and Kobayashi (2011).
///
/// Conventions, which differ between the two arguments and are the thing to get
/// right: residuals `u` arrive **k x tt**, one column per period, because that is
/// the shape `arma::reshape(y - z * a, k, tt)` produces. Latent scales `w` are
/// **tt x k**, one column per equation, because that is the orientation the
/// per-period precision is built from with `vectorise(trans(w))`. Neither is
/// arbitrary; both match what the samplers around them already hold.

namespace bayests::core
{

/// The two constants a quantile enters a model through.
struct AldShape
{
    double theta = 0.0; ///< The skew the mixture carries, (1 - 2q) / (q (1 - q)).
    double tau2 = 8.0;  ///< The variance multiplier, 2 / (q (1 - q)).
};

/// The shape constants of one quantile.
///
/// At q = 0.5 this is exactly {0, 8}: the skew vanishes and the asymmetric
/// Laplace becomes the symmetric one, which is the identity the unit test pins
/// and the case in which a model that has dropped the offset somewhere still
/// looks correct.
inline AldShape ald_shape(const double quantile)
{
    if (!(quantile > 0.0 && quantile < 1.0))
    {
        throw std::invalid_argument("the quantile must lie in (0, 1)");
    }

    const double denominator = quantile * (1.0 - quantile);
    AldShape shape;
    shape.theta = (1.0 - 2.0 * quantile) / denominator;
    shape.tau2 = 2.0 / denominator;
    return shape;
}

/// The log density of the asymmetric Laplace at residual `u` with scale `s`.
///
///     log f(u) = log q + log(1 - q) - log s - rho_q(u / s),
///     rho_q(v) = v (q - 1{v < 0}).
///
/// Closed form and marginal of the latent scale, which is why the log
/// likelihood of these models needs no drawn state and no determinant.
inline double ald_log_density(const double u, const double s, const double q)
{
    const double v = u / s;
    const double check = v * (q - (v < 0.0 ? 1.0 : 0.0));
    return std::log(q) + std::log(1.0 - q) - std::log(s) - check;
}

/// One draw of every latent scale, tt x k.
///
/// The full conditional is generalised inverse Gaussian at index 1/2. Writing
/// the kernel out,
///
///     p(w | u) ~ w^{-1/2} exp(-0.5 [ (u^2 / (tau2 s)) / w
///                                    + (theta^2 / (tau2 s) + 2 / s) w ]),
///
/// which is GIG(1/2, chi, psi) with chi = u^2 / (tau2 s) and
/// psi = (theta^2 + 2 tau2) / (tau2 s). A GIG at index +1/2 is the reciprocal of
/// one at -1/2, and GIG(-1/2, psi, chi) is the inverse Gaussian with mean
/// sqrt(psi / chi) and shape psi -- so the draw is one inverse Gaussian variate
/// inverted. Both parameters simplify:
///
///     mean  = sqrt(theta^2 + 2 tau2) / |u|,
///     shape = (theta^2 + 2 tau2) / (tau2 s).
///
/// @param u k x tt residuals, one column per period.
/// @param u_scale k scales, one per equation.
/// @param shape the quantile's constants.
///
/// A residual of exactly zero would send the mean to infinity, so `u` is floored
/// at a value far below any scale a model of standardised data works at. The
/// floor is on the magnitude only: the sign of a zero residual does not matter,
/// because the conditional depends on u through its square.
inline arma::mat draw_ald_weights(const arma::mat &u, const arma::vec &u_scale,
                                  const AldShape &shape)
{
    const arma::uword k = u.n_rows;
    const arma::uword tt = u.n_cols;

    const double numerator = shape.theta * shape.theta + 2.0 * shape.tau2;
    const double root_numerator = std::sqrt(numerator);

    // Floor chosen so that the resulting mean stays comfortably inside double
    // precision while being smaller than any residual a real model produces.
    constexpr double kResidualFloor = 1e-12;

    arma::vec mu(k * tt);
    arma::vec lambda(k * tt);
    for (arma::uword t = 0; t < tt; t++)
    {
        for (arma::uword i = 0; i < k; i++)
        {
            const double magnitude = std::max(std::abs(u(i, t)), kResidualFloor);
            const arma::uword pos = t * k + i;
            mu(pos) = root_numerator / magnitude;
            lambda(pos) = numerator / (shape.tau2 * u_scale(i));
        }
    }

    const arma::vec inverted = inverse_gaussian(mu, lambda);

    // Back to tt x k, which is the orientation the per-period precision is built
    // from.
    arma::mat w(tt, k);
    for (arma::uword t = 0; t < tt; t++)
    {
        for (arma::uword i = 0; i < k; i++)
        {
            w(t, i) = 1.0 / inverted(t * k + i);
        }
    }
    return w;
}

/// One draw of the scale of the asymmetric Laplace, k elements.
///
/// The residuals and the latent scales both carry information about it, so its
/// inverse gamma conditional collects two sums:
///
///     shape = shape_0 + 3 T / 2,
///     rate  = rate_0 + sum_t [ (u_t - theta w_t)^2 / (2 tau2 w_t) + w_t ].
///
/// The 3T/2 rather than T/2 is where the exponential prior on `w` enters -- it
/// contributes T of the T + T/2, and forgetting it biases the scale upwards.
///
/// @param u k x tt residuals, one column per period.
/// @param w tt x k latent scales, one column per equation.
/// @param post_shape the already updated shape, shape_0 + 3 T / 2.
/// @param prior_rate rate_0, one per equation.
inline arma::vec draw_ald_scale(const arma::mat &u, const arma::mat &w, const AldShape &shape,
                                const arma::vec &post_shape, const arma::vec &prior_rate)
{
    const arma::uword k = u.n_rows;
    const arma::uword tt = u.n_cols;

    arma::vec out(k);
    for (arma::uword i = 0; i < k; i++)
    {
        double rate = prior_rate(i);
        for (arma::uword t = 0; t < tt; t++)
        {
            const double weight = w(t, i);
            const double centred = u(i, t) - shape.theta * weight;
            rate += centred * centred / (2.0 * shape.tau2 * weight) + weight;
        }

        // Armadillo's randg takes a scale, so the rate is inverted at the point
        // of use -- the same arrangement every other gamma draw here makes.
        out(i) = 1.0 / arma::randg<double>(arma::distr_param(post_shape(i), 1.0 / rate));
    }
    return out;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_ALD_SUPPORT_H
