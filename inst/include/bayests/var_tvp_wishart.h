// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_TVP_WISHART_H
#define BAYESTS_VAR_TVP_WISHART_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR whose coefficients follow a random walk, with a Wishart prior on the
/// error precision and BVS variable selection.
///
/// The coefficient path is drawn as a block with the simulation smoother of
/// Durbin and Koopman (2002), so the whole path moves at once rather than one
/// period at a time.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
class VarTvpWishartSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarTvpWishartDraws draw_coefficients(const VarTvpWishartInput &input,
                                        Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, holding the
    /// coefficients and the precision at their last in-sample values:
    /// `draws.a` and `draws.u_sigma_inv` are expected to carry that period
    /// alone, one column per draw.
    ForecastDraws forecast(const VarTvpWishartInput &input,
                           const VarTvpWishartDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. `draws.a` carries the whole
    /// coefficient path -- every period is evaluated under its own
    /// coefficients -- while `draws.u_sigma_inv` carries the last period only,
    /// which is the precision this model scores every observation under.
    arma::mat log_likelihood(const VarTvpWishartInput &input,
                             const VarTvpWishartDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_TVP_WISHART_H
