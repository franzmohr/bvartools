// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_QUANTILE_GRID_H
#define BAYESTS_CORE_MODELS_QUANTILE_GRID_H

#include "bayests/arma.h"

#include <algorithm>
#include <cmath>

/// @file quantile_grid.h
/// @brief The distribution a grid of conditional quantiles describes.
///
/// A structural quantile VAR estimates each equation at a grid of quantile
/// levels tau_1 < ... < tau_Q. At one set of regressors that gives Q conditional
/// quantiles q_1, ..., q_Q, and this file turns them into a distribution the
/// model can be evaluated and simulated under:
///
///   - **Rearranged.** Estimated separately, the quantiles need not be
///     increasing in tau. Sorting them is the rearrangement of Chernozhukov,
///     Fernandez-Val and Galichon (2010), which can only move each estimate
///     closer to the true quantile function.
///   - **Linear between the grid points.** The quantile function is
///     interpolated linearly, so the density is constant on each interval,
///     (tau_{j+1} - tau_j) / (q_{j+1} - q_j).
///   - **Exponential beyond them.** Below q_1 the distribution function is
///     tau_1 exp(-lambda_L (q_1 - x)), above q_Q it is
///     1 - (1 - tau_Q) exp(-lambda_U (x - q_Q)). Each rate is the smaller --
///     the heavier tail -- of two: the one that makes the density continuous at
///     the end of the grid, which describes the data well where the outermost
///     quantiles are well estimated, and the rate of the asymmetric Laplace the
///     outermost level was estimated under, lambda_L = (1 - tau_1) / s_1 and
///     lambda_U = tau_Q / s_Q with s its scale. The second is estimated from
///     every residual of the sample and so stays bounded where the first
///     explodes: two outermost quantiles landing close together make the
///     continuous tail almost vertical, and a single observation beyond it would
///     otherwise score as all but impossible. Without scales the continuous
///     rate is used alone.
///
/// Neighbouring quantiles are kept at least kMinQuantileGap times the spread of
/// the grid apart, which keeps a density finite where two estimates coincide.
/// The two directions -- the log density of a value and the value at a level --
/// are exact inverses of each other, which is what the unit test pins.
///
/// Chavleishvili, S., & Manganelli, S. (2019). Forecasting and stress testing
/// with quantile vector autoregression. ECB Working Paper 2330.
///
/// Chernozhukov, V., Fernandez-Val, I., & Galichon, A. (2010). Quantile and
/// probability curves without crossing. Econometrica, 78(3), 1093-1125.

namespace bayests::core
{

constexpr double kMinQuantileGap = 1e-8;

/// The rearranged quantiles and the slopes of the distribution they describe.
struct QuantileGrid
{
    arma::vec tau;   ///< The levels, strictly increasing in (0, 1).
    arma::vec q;     ///< The quantiles at those levels, sorted.
    arma::vec gap;   ///< q_{j+1} - q_j, floored; Q - 1 elements.
    double lower = 0.0; ///< lambda_L, the rate of the lower tail.
    double upper = 0.0; ///< lambda_U, the rate of the upper tail.
};

/// Builds the distribution from one set of conditional quantiles, which need
/// not be sorted. `lower_scale` and `upper_scale` are the asymmetric Laplace
/// scales of the outermost levels, which bound the tail rates; zero leaves the
/// rates to continuity alone.
inline QuantileGrid quantile_grid(const arma::vec &tau, const arma::vec &quantiles,
                                  const double lower_scale = 0.0, const double upper_scale = 0.0)
{
    QuantileGrid g;
    g.tau = tau;
    g.q = arma::sort(quantiles);
    const arma::uword n = g.q.n_elem;
    // Coinciding quantiles are pushed apart by the floor, so that the density
    // and its inverse below see the same, strictly increasing, grid.
    const double floor = kMinQuantileGap * std::max(g.q(n - 1) - g.q(0), 1.0);
    for (arma::uword j = 1; j < n; j++)
    {
        g.q(j) = std::max(g.q(j), g.q(j - 1) + floor);
    }
    g.gap = arma::diff(g.q);

    // Continuity of the density at either end of the grid.
    const double first_density = (tau(1) - tau(0)) / g.gap(0);
    const double last_density = (tau(n - 1) - tau(n - 2)) / g.gap(n - 2);
    g.lower = first_density / tau(0);
    g.upper = last_density / (1.0 - tau(n - 1));
    if (lower_scale > 0.0 && upper_scale > 0.0)
    {
        g.lower = std::min(g.lower, (1.0 - tau(0)) / lower_scale);
        g.upper = std::min(g.upper, tau(n - 1) / upper_scale);
    }
    return g;
}

/// The log density at `x`.
inline double quantile_grid_log_density(const QuantileGrid &g, const double x)
{
    const arma::uword n = g.q.n_elem;
    if (x < g.q(0))
    {
        return std::log(g.tau(0) * g.lower) - g.lower * (g.q(0) - x);
    }
    if (x >= g.q(n - 1))
    {
        return std::log((1.0 - g.tau(n - 1)) * g.upper) - g.upper * (x - g.q(n - 1));
    }
    // The interval [q_j, q_{j+1}) holding x.
    const arma::uword j = static_cast<arma::uword>(
        std::upper_bound(g.q.begin(), g.q.end(), x) - g.q.begin()) - 1;
    return std::log((g.tau(j + 1) - g.tau(j)) / g.gap(j));
}

/// The value at level `u` in (0, 1): the inverse of the distribution function.
inline double quantile_grid_value(const QuantileGrid &g, const double u)
{
    const arma::uword n = g.q.n_elem;
    if (u < g.tau(0))
    {
        return g.q(0) + std::log(u / g.tau(0)) / g.lower;
    }
    if (u >= g.tau(n - 1))
    {
        return g.q(n - 1) - std::log((1.0 - u) / (1.0 - g.tau(n - 1))) / g.upper;
    }
    const arma::uword j = static_cast<arma::uword>(
        std::upper_bound(g.tau.begin(), g.tau.end(), u) - g.tau.begin()) - 1;
    return g.q(j) + (u - g.tau(j)) / (g.tau(j + 1) - g.tau(j)) * g.gap(j);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_QUANTILE_GRID_H
