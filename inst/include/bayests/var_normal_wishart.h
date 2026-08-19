// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_NORMAL_WISHART_H
#define BAYESTS_VAR_NORMAL_WISHART_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR with a normal prior on the coefficients and a Wishart prior on the
/// error precision, optionally with SSVS or BVS variable selection.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VarNormalWishartSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarNormalWishartDraws draw_coefficients(const VarNormalWishartInput &input,
                                            Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw.
    ///
    /// Requires `input.forecast.z` and, when the model has coefficients,
    /// `draws.a`. Throws std::invalid_argument if either is missing or if
    /// spec.h is zero.
    ForecastDraws forecast(const VarNormalWishartInput &input,
                           const VarNormalWishartDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior
    /// draw, one column per observation, as expected by WAIC and PSIS-LOO.
    arma::mat log_likelihood(const VarNormalWishartInput &input,
                             const VarNormalWishartDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_NORMAL_WISHART_H
