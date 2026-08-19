// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_NORMAL_WISHART_H
#define BAYESTS_VEC_NORMAL_WISHART_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC with a normal prior on the coefficients and a Wishart prior on the
/// error precision, optionally with SSVS or BVS variable selection.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecNormalWishartSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecNormalWishartDraws draw_coefficients(const VecNormalWishartInput &input,
                                            Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// The draws are rewritten in the level VAR parameterisation and the path is
    /// simulated by VarNormalWishartSampler, since a VEC and its level VAR are
    /// the same model and only one of them has a recursion worth writing twice.
    ///
    /// `input.forecast.z` is consequently expected in that level layout -- p + 1
    /// lags of the endogenous variables, s + 2 blocks of the unmodelled ones,
    /// then the unrestricted deterministic terms and the ones restricted to the
    /// cointegration space -- and not in the differenced layout `input.train.z`
    /// uses. The level history cannot be recovered from differenced regressors,
    /// so it has to be supplied rather than derived.
    ///
    /// Requires `input.forecast.z` and, when the model has coefficients,
    /// `draws.a`. Throws std::invalid_argument if either is missing, if spec.h
    /// is zero, or if the regressors do not match the converted coefficients.
    ForecastDraws forecast(const VecNormalWishartInput &input,
                           const VecNormalWishartDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior
    /// draw, one column per observation, as expected by WAIC and PSIS-LOO.
    arma::mat log_likelihood(const VecNormalWishartInput &input,
                             const VecNormalWishartDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_NORMAL_WISHART_H
