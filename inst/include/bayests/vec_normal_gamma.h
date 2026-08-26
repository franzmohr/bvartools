// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_NORMAL_GAMMA_H
#define BAYESTS_VEC_NORMAL_GAMMA_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC with a normal prior on the coefficients, the cointegration space prior of
/// Koop, Leon-Gonzalez and Strachan (2010) on beta, and independent gamma priors
/// on the error precisions, optionally with a constant covariance block and
/// SSVS or BVS variable selection.
///
/// VecNormalWishart's two coefficient blocks against VarNormalGamma's error
/// block. The only place the swap shows is Sigma's own posterior: the
/// cointegration space prior conditions alpha on Sigma, so the Wishart absorbs
/// a term back into its scale and its degrees of freedom, while independent
/// gammas have no conjugate update for it. bvartools' .bvecalg does not attempt
/// one either, and this follows it -- see the comment at the error block.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecNormalGammaSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecNormalGammaDraws draw_coefficients(const VecNormalGammaInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// The draws are rewritten in the level VAR parameterisation and the path is
    /// simulated by VarNormalWishartSampler, for the reason
    /// VecNormalWishartSampler::forecast() gives: a VEC and its level VAR are
    /// the same model, and only one of them has a recursion worth writing twice.
    /// `input.forecast.z` is consequently expected in that level layout, not in
    /// the differenced one `input.train.z` uses.
    ForecastDraws forecast(const VecNormalGammaInput &input,
                           const VecNormalGammaDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    arma::mat log_likelihood(const VecNormalGammaInput &input,
                             const VecNormalGammaDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_NORMAL_GAMMA_H
