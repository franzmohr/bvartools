// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_NORMAL_GAMMA_H
#define BAYESTS_VAR_NORMAL_GAMMA_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR with a normal prior on the coefficients and independent gamma priors on
/// the error precisions, optionally with a constant covariance block and with
/// SSVS or BVS variable selection.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VarNormalGammaSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarNormalGammaDraws draw_coefficients(const VarNormalGammaInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw.
    ///
    /// Uses `draws.a` only when `input.forecast.z` is present or the model is
    /// structural; without either, the path is the error process alone. In a
    /// structural model the last k(k-1)/2 coefficients are read as the
    /// contemporaneous matrix and applied to each period's draw.
    ForecastDraws forecast(const VarNormalGammaInput &input,
                           const VarNormalGammaDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior
    /// draw, one column per observation, as expected by WAIC and PSIS-LOO.
    arma::mat log_likelihood(const VarNormalGammaInput &input,
                             const VarNormalGammaDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_NORMAL_GAMMA_H
