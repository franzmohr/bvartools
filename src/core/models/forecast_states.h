// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_FORECAST_STATES_H
#define BAYESTS_CORE_MODELS_FORECAST_STATES_H

#include "bayests/arma.h"
#include "bayests/spec.h"

#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Carrying a random walk over a forecast horizon.
///
/// A time-varying model's coefficients, covariance block and log-volatilities
/// are random walks, so the model's own forecast of any of them at step i is
/// the value at the end of the sample plus i innovations. Holding them at that
/// value -- ForecastStates::hold -- forecasts from a model in which the drift
/// stops where the sample does. What that leaves out is more than width: the
/// variance a log-volatility implies grows as exp(i sigma / 2) in expectation,
/// so a held volatility understates it on average, and a product of drifting
/// coefficient matrices has a different mean from the product of held ones.
///
/// The state at step i is the one the observation at step i is generated from,
/// so the first innovation is added before the first horizon, not after it:
/// y_{T+1} is drawn under a_{T+1} = a_T + eta, exactly as y_t is under a_t in
/// sample.

inline bool simulates_states(const VarSpec &spec)
{
    return spec.forecast_states == ForecastStates::simulate;
}

/// Refuses draws of a random walk's innovation variances that a forecast cannot
/// simulate from: one row per element of the state, one column per draw, and
/// nothing negative, which would reach the square root below as a NaN.
inline void require_state_variances(const arma::mat &sigma, const arma::uword n_state,
                                    const arma::uword draws, const std::string &name)
{
    if (sigma.n_elem == 0)
    {
        throw std::invalid_argument(
            "simulating " + name + " forward over the forecast horizon needs posterior draws of "
            "its innovation variances, which are missing; forecast_states = hold forecasts "
            "without them");
    }
    if (sigma.n_rows != n_state || sigma.n_cols != draws)
    {
        throw std::invalid_argument(
            "posterior draws of the innovation variances of " + name + " must be " +
            std::to_string(n_state) + " x " + std::to_string(draws) + ", got " +
            std::to_string(sigma.n_rows) + " x " + std::to_string(sigma.n_cols));
    }
    if (sigma.min() < 0.0)
    {
        throw std::invalid_argument("posterior draws of the innovation variances of " + name +
                                    " hold a negative value");
    }
}

/// Refuses the value a state starts the horizon from unless it is one period of
/// `rows` per draw.
inline void require_period_draws(const arma::mat &values, const arma::uword rows,
                                 const arma::uword draws, const std::string &name)
{
    if (values.n_rows != rows || values.n_cols != draws)
    {
        throw std::invalid_argument(
            "simulating the states forward over the forecast horizon needs posterior draws of " +
            name + " at the last in-sample period, " + std::to_string(rows) +
            " per draw over " + std::to_string(draws) + " draws, got " +
            std::to_string(values.n_rows) + " x " + std::to_string(values.n_cols));
    }
}

/// Refuses a selection mask of the wrong height. Empty is no selection.
inline void require_state_mask(const arma::mat &lambda, const arma::uword rows,
                               const arma::uword draws, const std::string &name)
{
    if (lambda.n_elem > 0 && (lambda.n_rows != rows || lambda.n_cols != draws))
    {
        throw std::invalid_argument(
            "posterior draws of the inclusion indicators of " + name + " must be " +
            std::to_string(rows) + " x " + std::to_string(draws) + ", got " +
            std::to_string(lambda.n_rows) + " x " + std::to_string(lambda.n_cols));
    }
}

/// One step of a random walk: `state += sqrt(sigma) % N(0, I)`.
///
/// A position `mask` switches off gets no innovation. BVS stores a time-varying
/// path with an excluded regressor's row zeroed in every period, and the model
/// is y = Z Lambda theta: an excluded coefficient that drifted away from zero
/// over the horizon would put back a regressor the posterior draw left out.
/// Standard normals are drawn for every position either way, so a mask does not
/// change which draws the rest of the forecast gets.
inline void step_random_walk(arma::vec &state, const arma::vec &sigma, const arma::vec &mask)
{
    arma::vec innovation = arma::sqrt(sigma) % arma::randn<arma::vec>(state.n_elem);
    if (mask.n_elem > 0)
    {
        innovation %= mask;
    }
    state += innovation;
}

/// The strict lower triangle of `matrix`, row by row -- the inverse of
/// fill_strict_lower_triangle(), and so the order Psi's free elements, their
/// innovation variances and their inclusion indicators are all counted in.
/// Used to recover the packed form from a posterior that stores Psi whole.
inline arma::vec pack_strict_lower_triangle(const arma::mat &matrix)
{
    const arma::uword k = matrix.n_rows;
    arma::vec packed(k * (k - 1) / 2);
    for (arma::uword i = 1; i < k; i++)
    {
        packed.subvec(i * (i - 1) / 2, (i + 1) * i / 2 - 1) =
            arma::trans(matrix.submat(i, 0, i, i - 1));
    }
    return packed;
}

/// The symmetric square root of the covariance whose precision is
/// `precision`, which is what a forecast error is drawn through: root * N(0, I).
/// The same factorisation every forecast here uses, for a precision that is
/// rebuilt at each horizon because something in it drifts.
inline arma::mat covariance_root(const arma::mat &precision)
{
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec,
                  arma::solve(precision, arma::eye<arma::mat>(precision.n_rows, precision.n_rows)));
    return eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::trans(eigvec);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_FORECAST_STATES_H
