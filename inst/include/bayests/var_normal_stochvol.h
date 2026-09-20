// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_NORMAL_STOCHVOL_H
#define BAYESTS_VAR_NORMAL_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR with a normal prior on the coefficients and stochastic volatility in
/// the errors, optionally with a covariance block and BVS variable selection.
///
/// The log-volatility is drawn with the ten-component normal mixture of Omori
/// et al. (2007), which turns the non-linear measurement equation into a
/// conditionally linear one.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VarNormalStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarNormalStochvolDraws draw_coefficients(const VarNormalStochvolInput &input,
                                             Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw from the last in-sample
    /// volatility: `draws.u_sigma_inv` (k * k) and `draws.u_omega_inv` (k) are
    /// expected to carry that period alone, one column per draw.
    ///
    /// Under ForecastStates::simulate, the default, the log-volatilities take
    /// one step of their random walk per horizon, by `draws.h_sigma`, and the
    /// precision is rebuilt from them and the constant `draws.psi` at every
    /// horizon. Under ForecastStates::hold the volatility stays where the
    /// sample ends and only `draws.u_sigma_inv` is read of the two.
    ForecastDraws forecast(const VarNormalStochvolInput &input,
                           const VarNormalStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. Each period is evaluated
    /// under its own precision matrix, so `draws.u_sigma_inv` carries the full
    /// (k * k * tt) path here rather than the single period the forecast uses.
    arma::mat log_likelihood(const VarNormalStochvolInput &input,
                             const VarNormalStochvolDraws &draws) const;

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
    arma::mat predictive_log_density(const VarNormalStochvolInput &input,
                                     const VarNormalStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_NORMAL_STOCHVOL_H
