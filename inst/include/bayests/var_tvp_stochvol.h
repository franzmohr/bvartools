// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_TVP_STOCHVOL_H
#define BAYESTS_VAR_TVP_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR whose coefficients follow a random walk, with stochastic volatility in
/// the errors, an optional time-varying covariance block and BVS variable
/// selection.
///
/// This is the two moving parts of VarTvpGamma and VarNormalStochvol in one
/// model: the coefficient path is drawn as a block with the simulation smoother
/// of Durbin and Koopman (2002), and the log-volatility with the ten-component
/// normal mixture of Omori et al. (2007).
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VarTvpStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarTvpStochvolDraws draw_coefficients(const VarTvpStochvolInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw from the last in-sample
    /// period: `draws.a`, `draws.u_sigma_inv`, `draws.u_omega_inv` (k) and
    /// `draws.psi` (k * k) are expected to carry that period alone, one column
    /// per draw.
    ///
    /// Under ForecastStates::simulate, the default, the coefficients, Psi and
    /// the log-volatilities each take one step of their random walk per
    /// horizon, by `draws.a_sigma`, `draws.psi_sigma` and `draws.h_sigma`, with
    /// the positions `draws.a_lambda` and `draws.psi_lambda` exclude held at
    /// zero, and the precision is rebuilt from them at every horizon. Under
    /// ForecastStates::hold all of it stays where the sample ends and only
    /// `draws.a` and `draws.u_sigma_inv` are read.
    ForecastDraws forecast(const VarTvpStochvolInput &input,
                           const VarTvpStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. `draws.a` carries the whole
    /// coefficient path and `draws.u_sigma_inv` the whole precision path, one
    /// k x k block per period: every period is evaluated under its own
    /// coefficients and its own precision.
    arma::mat log_likelihood(const VarTvpStochvolInput &input,
                             const VarTvpStochvolDraws &draws) const;

    /// The log predictive density of what the horizon realised, draws x scored
    /// periods -- one row per posterior draw, one column per period of
    /// `input.test.y`.
    ///
    /// Each column conditions on the realised observations before it rather than
    /// on a simulated path, so the log of the mean of exp() over draws is the
    /// one step ahead predictive density given everything known up to that
    /// period, and those sum over the horizons to the log predictive likelihood
    /// of the whole realised path.
    ///
    /// Expects the draws a forecast expects -- the last in-sample period of
    /// what moves, plus the innovation variances to move it by -- and carries
    /// them forward one step per scored period, honouring
    /// `spec.forecast_states`. Under `simulate` that is one sampled state path
    /// per draw, so the score is drawn rather than computed and `/model/seed`
    /// is what repeats it. Requires `input.test.y` and `input.forecast.x`;
    /// throws for a structural model, which this expression is not the density
    /// of.
    arma::mat predictive_log_density(const VarTvpStochvolInput &input,
                                     const VarTvpStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_TVP_STOCHVOL_H
