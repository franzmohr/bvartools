// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_TVP_STOCHVOL_H
#define BAYESTS_VEC_TVP_STOCHVOL_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC whose coefficients follow a random walk -- the cointegration vectors
/// among them -- with stochastic volatility in the errors, an optional
/// time-varying covariance block and BVS variable selection.
///
/// This is VarTvpStochvol plus one Gibbs block, and the block is the whole
/// difference. A VEC's first n_alpha regressors are beta' w_{t-1}, so they are
/// not data: they are rebuilt from the current draw of beta. The sampler
/// therefore alternates between two state paths conditioned on each other --
/// `a` given the regressors beta implies, then beta given the loadings `a`
/// carries -- each drawn as a block with the simulation smoother of Durbin and
/// Koopman (2002), with the log-volatility drawn from the ten-component normal
/// mixture of Omori et al. (2007) as in the VAR.
///
/// The alternation, and the shape of both blocks, follows bvartools'
/// .bvectvpalg. Two things it does are settled differently here, and both are
/// documented where they happen: the coefficient blocks are drawn with the
/// smoother rather than as one (n_a tt)-dimensional normal, and BVS may not
/// reach the loadings.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecTvpStochvolSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecTvpStochvolDraws draw_coefficients(const VecTvpStochvolInput &input,
                                          Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// Everything this model estimates moves with time, so the forecast starts
    /// from the last in-sample period of all of it: `draws.a`, `draws.beta` and
    /// `draws.u_sigma_inv` are expected to carry that period alone, one column
    /// per draw. Those are then rewritten in the level VAR parameterisation and
    /// the path is simulated by VarNormalWishartSampler, for the reason
    /// VecNormalWishartSampler::forecast() gives: a VEC and its level VAR are
    /// the same model, and only one of them has a recursion worth writing twice.
    ///
    /// `input.forecast.z` is consequently expected in that level layout -- p + 1
    /// lags of the endogenous variables, s + 2 blocks of the unmodelled ones,
    /// then the unrestricted deterministic terms and the ones restricted to the
    /// cointegration space -- and not in the differenced layout `input.train.z`
    /// uses. The level history cannot be recovered from differenced regressors,
    /// so it has to be supplied rather than derived.
    ForecastDraws forecast(const VecTvpStochvolInput &input,
                           const VecTvpStochvolDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. `draws.a` and `draws.beta`
    /// carry the whole path -- every period is evaluated under its own
    /// coefficients and its own cointegration vectors -- while
    /// `draws.u_sigma_inv` carries the last period only, which is the precision
    /// this model scores every observation under.
    arma::mat log_likelihood(const VecTvpStochvolInput &input,
                             const VecTvpStochvolDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_TVP_STOCHVOL_H
