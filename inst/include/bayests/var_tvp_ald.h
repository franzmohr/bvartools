// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_TVP_ALD_H
#define BAYESTS_VAR_TVP_ALD_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR estimated at a conditional quantile whose coefficients follow a random
/// walk, with optional BVS variable selection.
///
/// VarNormalAld with the coefficients turned into a path, exactly as
/// VarTvpStochvol is VarNormalStochvol with the same change, and every word of
/// VarNormalAld's class comment still applies -- the asymmetric Laplace scale
/// mixture, the absent covariance block, the absent forecast, and the working
/// likelihood the intervals should not be read as a posterior over.
///
/// What the drift buys is the thing a constant-coefficient quantile model cannot
/// say. A quantile regression already lets the spread of the errors differ from
/// the middle of the distribution; letting the coefficients move as well is what
/// separates "the tail became noisier" from "the relationship at the tail
/// changed" -- and those are different claims about the same widening fan.
///
/// Seven blocks against five. The coefficient path is drawn whole with the
/// simulation smoother of Durbin and Koopman (2002), against the per-period
/// covariance tau2 * s_i * w_it and the offset theta * w_it; then its state
/// innovation variance and its pre-sample state; then the latent scales, the
/// asymmetric Laplace scale, and the rebuild.
///
/// The measurement offset is the one thing easy to lose here, because it has to
/// reach three places that each look self-contained: the response handed to the
/// smoother, the residual the BVS sweep scores, and the log likelihood. A model
/// that forgets it in any one of them is estimating the median wherever it
/// forgot, and at q = 0.5 -- where theta is zero -- nothing would show.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG.
///
/// Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation
/// smoother for state space time series analysis. Biometrika, 89(3), 603-615.
///
/// Kozumi, H., & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
/// quantile regression. Journal of Statistical Computation and Simulation,
/// 81(11), 1565-1578.
class VarTvpAldSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarTvpAldDraws draw_coefficients(const VarTvpAldInput &input, Reporter &reporter) const;

    /// Always throws. A quantile VAR has no forecast: see VarNormalAldSampler.
    ForecastDraws forecast(const VarTvpAldInput &input, const VarTvpAldDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods. The asymmetric Laplace density
    /// under each period's own coefficient vector.
    arma::mat log_likelihood(const VarTvpAldInput &input, const VarTvpAldDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_TVP_ALD_H
