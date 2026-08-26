// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_TVP_WISHART_H
#define BAYESTS_VEC_TVP_WISHART_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC whose coefficients follow a random walk -- the cointegration vectors
/// among them -- with a Wishart prior on the error precision and BVS variable
/// selection. It carries no psi block: the error covariance is the Wishart
/// precision alone.
///
/// VecTvpStochvol's two state blocks against VarTvpWishart's error block. The
/// alternation is the one VecTvpStochvol describes -- `a` given the regressors
/// beta implies, then beta given the loadings `a` carries, each drawn with the
/// simulation smoother of Durbin and Koopman (2002) -- and only what the two are
/// conditioned on differs, being the same k x k precision in every period here.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecTvpWishartSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecTvpWishartDraws draw_coefficients(const VecTvpWishartInput &input,
                                         Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// `draws.a` and `draws.beta` are expected to carry the last in-sample
    /// period alone, one column per draw; the precision does not move, so it is
    /// passed whole. The draws are then rewritten in the level VAR
    /// parameterisation and simulated by VarNormalWishartSampler, so
    /// `input.forecast.z` is expected in that level layout and not in the
    /// differenced one `input.train.z` uses.
    ForecastDraws forecast(const VecTvpWishartInput &input,
                           const VecTvpWishartDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. `draws.a` and `draws.beta`
    /// carry the whole path -- every period is evaluated under its own
    /// coefficients and its own cointegration vectors -- against the single
    /// precision of each draw.
    arma::mat log_likelihood(const VecTvpWishartInput &input,
                             const VecTvpWishartDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_TVP_WISHART_H
