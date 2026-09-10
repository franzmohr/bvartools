// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_NORMAL_ALD_H
#define BAYESTS_VAR_NORMAL_ALD_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR estimated at a conditional quantile, with a normal prior on the
/// coefficients and optional BVS variable selection.
///
/// Minimising the quantile loss at q is maximising the likelihood of an
/// asymmetric Laplace distribution, and that distribution is a scale mixture of
/// normals. With
///
///     theta = (1 - 2q) / (q (1 - q)),   tau2 = 2 / (q (1 - q)),
///
/// the model
///
///     y_it = x_t' a_i + theta w_it + e_it,  e_it ~ N(0, tau2 s_i w_it),
///     w_it ~ Exp(1 / s_i),
///
/// has the q-th conditional quantile of y_it at x_t' a_i. That representation
/// is what makes this a Gibbs sampler rather than an optimiser: conditional on
/// the latent scales w, every block is one this library already draws.
///
/// Five blocks. The coefficients, from the same weighted normal posterior
/// VarNormalStochvol draws -- only the per-period variance comes from
/// tau2 * s_i * w_it rather than from a log-volatility, and the response
/// carries the offset theta * w_it. The latent scales, whose full conditional
/// is generalised inverse Gaussian at index 1/2, drawn as the reciprocal of an
/// inverse Gaussian. The scale s_i, from an inverse gamma. Then the precision
/// is rebuilt and, under selection, the BVS sweep runs against the same
/// offset residual.
///
/// Three things this model does not have, each for a reason worth knowing.
///
///   - **No covariance block.** Psi is a triangular rotation of the errors.
///     Conditional on w the equations are independent, and rotating them is
///     exactly what stops the estimand being a quantile: the rotated residual
///     is a combination of equations, and the q-th quantile of a combination is
///     not the combination of q-th quantiles. validate() rejects `covar`.
///   - **No forecast.** The h step quantile is not the quantile of the iterated
///     one step quantiles, so there is no number to return that could be read as
///     one. validate() rejects a non-zero horizon rather than letting the
///     front-end produce something that looks like a forecast, and forecast()
///     throws if it is reached anyway.
///   - **No calibrated intervals.** The asymmetric Laplace is a working
///     likelihood, not a claim about the data. The posterior locates the
///     quantile, but its spread needs the sandwich adjustment of Yang, Wang and
///     He (2016), which this does not apply. Read the spread as a diagnostic,
///     not as a credible interval.
///
/// Structural models *are* allowed: contemporaneous terms are regressors like
/// any other, and adding them to an equation leaves its quantile reading intact.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG.
///
/// Kozumi, H., & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian
/// quantile regression. Journal of Statistical Computation and Simulation,
/// 81(11), 1565-1578.
///
/// Yang, Y., Wang, H. J., & He, X. (2016). Posterior inference in Bayesian
/// quantile regression with asymmetric Laplace likelihood. International
/// Statistical Review, 84(3), 327-344.
class VarNormalAldSampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarNormalAldDraws draw_coefficients(const VarNormalAldInput &input,
                                        Reporter &reporter) const;

    /// Always throws. A quantile VAR has no forecast: see the class comment.
    ///
    /// Declared so that the front-ends and any embedded host can call the same
    /// three methods on every sampler, and so that the refusal carries a reason
    /// rather than being a missing symbol.
    ForecastDraws forecast(const VarNormalAldInput &input, const VarNormalAldDraws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior draw,
    /// one column per observation, as expected by WAIC and PSIS-LOO.
    ///
    /// The asymmetric Laplace density itself, which is closed form and marginal
    /// of the latent scales. That is what makes this simpler than the stochastic
    /// volatility models' version rather than harder: there is no determinant to
    /// recompute per period, and no conditioning on a drawn state.
    arma::mat log_likelihood(const VarNormalAldInput &input,
                             const VarNormalAldDraws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_NORMAL_ALD_H
