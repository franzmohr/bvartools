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

    /// Simulates one forecast path per posterior draw, holding the volatility
    /// at its last in-sample value: `draws.u_sigma_inv` is expected to carry
    /// that period alone, one vectorised k x k precision per column.
    ForecastDraws forecast(const VarNormalStochvolInput &input,
                           const VarNormalStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. Each period is evaluated
    /// under its own precision matrix, so `draws.u_sigma_inv` carries the full
    /// (k * k * tt) path here rather than the single period the forecast uses.
    arma::mat log_likelihood(const VarNormalStochvolInput &input,
                             const VarNormalStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_NORMAL_STOCHVOL_H
