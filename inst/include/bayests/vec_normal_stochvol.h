// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_NORMAL_STOCHVOL_H
#define BAYESTS_VEC_NORMAL_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests

{

/// VEC with a normal prior on the coefficients, the cointegration space prior of
/// Koop, Leon-Gonzalez and Strachan (2010) on beta, and stochastic volatility in
/// the errors, optionally with a constant covariance block and BVS variable
/// selection.
///
/// VecNormalWishart's two coefficient blocks against VarNormalStochvol's error
/// block, and the combination costs one thing beyond the substitution. The
/// cointegration space prior conditions alpha on the error precision, and here
/// there is a different one in every period; the prior takes their average over
/// the sample, which is what bvartools' .bvecalg does with its `g_i`. Beta's own
/// posterior then uses the per-period precisions in full, so the average appears
/// only where the prior does.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecNormalStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecNormalStochvolDraws draw_coefficients(const VecNormalStochvolInput &input,
                                             Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// The coefficients are constant, so only the precision has a period to be
    /// held at: `draws.u_sigma_inv` is expected to carry the last in-sample one
    /// alone. The draws are then rewritten in the level VAR parameterisation and
    /// simulated by VarNormalWishartSampler, so `input.forecast.z` is expected in
    /// that level layout and not in the differenced one `input.train.z` uses.
    ForecastDraws forecast(const VecNormalStochvolInput &input,
                           const VecNormalStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. Every period is scored under
    /// its own precision and under the one constant set of coefficients.
    arma::mat log_likelihood(const VecNormalStochvolInput &input,
                             const VecNormalStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_NORMAL_STOCHVOL_H
