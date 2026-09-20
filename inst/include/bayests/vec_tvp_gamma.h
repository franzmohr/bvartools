// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_TVP_GAMMA_H
#define BAYESTS_VEC_TVP_GAMMA_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC whose coefficients follow a random walk -- the cointegration vectors
/// among them -- with independent gamma priors on the error precisions, an
/// optional time-varying covariance block and BVS variable selection.
///
/// VecTvpStochvol's two state blocks against VarTvpGamma's error block. The
/// alternation is the one VecTvpStochvol describes -- `a` given the regressors
/// beta implies, then beta given the loadings `a` carries, each drawn with the
/// simulation smoother of Durbin and Koopman (2002).
///
/// The precision is constant across periods unless the model has a covariance
/// block, in which case Psi is a path and Psi_t' Omega^-1 Psi_t is one too. That
/// is the only reason this model's output shape depends on `covar`, and it is
/// the same dependence VarTvpGamma has.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VecTvpGammaSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecTvpGammaDraws draw_coefficients(const VecTvpGammaInput &input,
                                       Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// `draws.a` and `draws.beta` are expected to carry the last in-sample
    /// period alone, one column per draw, and so are `draws.u_sigma_inv` and
    /// `draws.psi` (k * k) when the model has a covariance block to make them
    /// move. `input.forecast.x` is expected in the level layout and not in the
    /// differenced one `input.train.z` uses.
    ///
    /// Under ForecastStates::simulate, the default, each horizon the
    /// coefficients step by `draws.a_sigma`, the cointegration vectors by their
    /// state equation as in VecTvpWishartSampler::forecast(), and Psi by
    /// `draws.psi_sigma`, with the positions the two selection masks exclude
    /// held at zero; the level VAR and the precision, rebuilt around the constant
    /// `draws.u_omega_inv` (k), are formed again at every horizon. Under
    /// ForecastStates::hold the draws are rewritten in the level VAR
    /// parameterisation once and simulated by VarNormalWishartSampler.
    ForecastDraws forecast(const VecTvpGammaInput &input,
                           const VecTvpGammaDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. `draws.a` and `draws.beta`
    /// carry the whole path; `draws.u_sigma_inv` carries one k x k block per
    /// period when the precision moves and the single matrix of each draw when
    /// it does not, and every period is scored under its own.
    arma::mat log_likelihood(const VecTvpGammaInput &input,
                             const VecTvpGammaDraws &draws) const;

    /// The log predictive density of what the horizon realised, draws x scored
    /// periods -- one row per posterior draw, one column per period of
    /// `input.test.y`, which for a VEC holds the realised levels.
    ///
    /// Scored in the level parameterisation this model forecasts in, under the
    /// states each draw reaches by stepping forward, with the change of basis
    /// made again at every period because the level coefficients are not linear
    /// in the states. Under `hold` nothing steps and the conversion is made
    /// once. See VarNormalWishartSampler::predictive_log_density() for what the
    /// number is and what it conditions on.
    arma::mat predictive_log_density(const VecTvpGammaInput &input,
                                     const VecTvpGammaDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_TVP_GAMMA_H
