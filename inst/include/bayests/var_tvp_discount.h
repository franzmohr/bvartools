// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VAR_TVP_DISCOUNT_H
#define BAYESTS_VAR_TVP_DISCOUNT_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VAR whose coefficients follow a random walk and whose error covariance
/// drifts, estimated in closed form rather than sampled.
///
/// The model is the matrix normal dynamic linear model of West and Harrison
/// (1997, ch. 16) with Uhlig's (1997) discounted Wishart on the error
/// precision. Two things make the posterior conjugate where a general
/// time varying parameter VAR's is not:
///
/// - Every equation of a VAR has the same regressors, so the coefficient
///   covariance factorises as `C_t` kron `Sigma_t` and only the n_reg square
///   factor is carried. That is also why the forecast scale is a scalar and no
///   (n_reg k) square covariance is ever formed.
/// - Neither the drift nor the volatility has a free variance parameter. The
///   drift is the discount `VarSpec::delta_beta`, which is a model rather than
///   a plug in: `R_t = C_{t-1} / delta` is exactly the predicted covariance of
///   a random walk whose innovation covariance is
///   `((1 - delta) / delta) C_{t-1}`, and that depends on the data only through
///   the past. The volatility is `VarSpec::delta_sigma` discounting the Wishart
///   posterior each period, which is a stochastic volatility law in its own
///   right -- multiplicative beta shocks to the precision, where
///   `stochvol_ocsn_2007` approximates a Gaussian autoregression in the log
///   variance.
///
/// **This is the one model in the project that is not a sampler**, and it is
/// the reason `VarTvpDiscountPosterior` holds a posterior rather than draws.
/// Setting both discount factors to one reproduces the conjugate normal
/// inverse Wishart posterior of a constant coefficient VAR exactly, which is
/// the check the unit test pins it against. Nothing here consumes the RNG, so
/// two runs on the same input agree to the bit and `/model/seed` has nothing to
/// repeat.
///
/// What it gives up against the samplers beside it is stated rather than
/// hidden. The volatility does not mean revert, and one discount factor stands
/// in for both the persistence and the variance of the volatility process, so
/// `delta_sigma` is not the autoregressive parameter of a log volatility. The
/// degrees of freedom converge to `1 / (1 - delta_sigma)` whatever the prior
/// says, so a short memory cannot carry a wide system.
///
/// Values in, values out: no files, no console, no global state. The reporter
/// is taken for cancellation and for one progress step per period, not per
/// draw -- there are none.
class VarTvpDiscountEstimator
{
public:
    /// Runs the filter and the retrospective pass. Reports progress once per
    /// period and honours an interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VarTvpDiscountPosterior estimate(const VarTvpDiscountInput &input,
                                     Reporter &reporter) const;

    /// I.i.d. draws from the posterior of one period, `draws` of them: first
    /// `Sigma ~ IW(df, df S_t)`, then `vec(Theta) | Sigma ~ N(m_t, Sigma kron
    /// C_t)`. Returns nparams x draws, the layout every sampler in this project
    /// returns, so a host can hand the result to the same consumers.
    ///
    /// Correct for that period alone. The smoothed posterior is not
    /// independent across periods, so calling this per period and joining the
    /// columns does not give a draw of the coefficient path.
    arma::mat draw_period(const VarTvpDiscountPosterior &posterior,
                          arma::uword period, arma::uword draws) const;

    /// Simulates one forecast path per draw from the last in-sample period,
    /// the coefficients taking one step of their random walk per horizon under
    /// the discount. `draws` decides how many, since the posterior it starts
    /// from carries none.
    ForecastDraws forecast(const VarTvpDiscountInput &input,
                           const VarTvpDiscountPosterior &posterior,
                           arma::uword draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, 1 x periods rather than draws x periods: the
    /// parameters are integrated out exactly, so there is one number per period
    /// and no Monte Carlo average to take over.
    ///
    /// A consumer that wants the draws x periods matrix WAIC and PSIS-LOO
    /// expect builds it with `draw_period()` instead, and pays for the draws it
    /// then averages over.
    arma::mat log_likelihood(const VarTvpDiscountInput &input,
                             const VarTvpDiscountPosterior &posterior) const;

    /// The log predictive density of what the horizon realised, 1 x scored
    /// periods. Requires `input.test.y` and `input.forecast.x`.
    ///
    /// Each period conditions on the realised observations before it, which
    /// this model gets for nothing: a filter is recursive, so the one step
    /// ahead predictive of every in-sample period is already in
    /// `VarTvpDiscountPosterior::loglik` and an expanding window needs no
    /// re-estimation at all. This entry point is for periods past the sample.
    arma::mat predictive_log_density(const VarTvpDiscountInput &input,
                                     const VarTvpDiscountPosterior &posterior) const;
};

} // namespace bayests

#endif // BAYESTS_VAR_TVP_DISCOUNT_H
