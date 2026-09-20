// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_TVP_DISCOUNT_H
#define BAYESTS_VEC_TVP_DISCOUNT_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// VEC whose loadings and short-run coefficients follow a random walk and whose
/// error covariance drifts, over a cointegration space that is **held fixed**,
/// estimated in closed form rather than sampled.
///
/// VarTvpDiscountEstimator over a design with the error correction term in
/// front of it. Everything that class says about the recursion applies here
/// unchanged, because conditional on `beta` a VEC *is* a VAR in the regressors
///
///     x_t = [ beta' w_t , the lagged differences, the deterministic terms ],
///
/// which every equation shares. That is the whole of the model, and the whole
/// of what it costs.
///
/// **Why the space cannot move here**, when the three sampling VECs beside it
/// let it. The measurement of a VEC is `alpha_t beta_t' w_t`, two latent blocks
/// multiplying each other, and a dynamic linear model needs its design known
/// before the period it explains. Conditioning on one of the two restores that,
/// which is why VecTvpWishartSampler and its siblings alternate two simulation
/// smoother passes -- and why a discount cannot replace the alternation, only
/// one pass inside it. Three things stand in the way of discounting the space
/// itself, and the third is decisive:
///
/// - Given the loadings, the design for vec(beta_t) is `kron(alpha, w_t')`,
///   which mixes the equations. The Kronecker factorisation that keeps this
///   model to an n_design square factor is gone, and the state covariance is a
///   full (k_beta rank) square object.
/// - The posterior would no longer be conjugate jointly with the discounted
///   Wishart on Sigma, so the closed form -- and with it the exact marginal
///   likelihood, the absence of a seed, and the absence of a burn-in -- goes
///   with it.
/// - TvpCointSpacePrior fixes the innovation variance of beta_t at the identity
///   and says why: only `alpha beta'` is identified, and the unit variance is
///   what pins beta's scale so that the level of the relation lives in alpha. A
///   discount replaces a fixed innovation variance with
///   `((1 - delta) / delta) C_{t-1}`, which moves with the data and with the
///   current uncertainty. That unpins the scale in the one place nothing else
///   pins it, and leaves a model that runs and produces plausible numbers while
///   alpha and beta slide along a ridge against each other.
///
/// So this is not a cheaper Koop, Leon-Gonzalez and Strachan (2011); it answers
/// the narrower question of how the adjustment to a given long-run relation has
/// moved. Whether the relation itself moved is what those three samplers are
/// for.
///
/// What the fixed space buys back is the marginal likelihood. The sum of
/// VecTvpDiscountPosterior::loglik is `p(Y | beta, rank, delta_beta,
/// delta_sigma)` exactly, at the cost of one pass over the sample, so a host can
/// compare candidate spaces, ranks and discounts against each other without
/// running a chain for any of them.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG, which only the forecast and the per-period draws touch --
/// estimate() consumes no random numbers at all.
class VecTvpDiscountEstimator
{
public:
    /// Runs the filter and the retrospective pass. Reports progress once per
    /// period and honours an interrupt thrown from the reporter.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecTvpDiscountPosterior estimate(const VecTvpDiscountInput &input, Reporter &reporter) const;

    /// I.i.d. draws from the posterior of one period, `draws` of them: first
    /// `Sigma ~ IW(df, df S_t)`, then `vec(Theta) | Sigma ~ N(m_t, Sigma kron
    /// C_t)`. Returns nparams x draws in the order a VEC's `a` is stored in,
    /// vec(alpha) first, so the result goes to the same consumers a sampled
    /// VEC posterior does.
    ///
    /// Correct for that period alone; see VarTvpDiscountEstimator::draw_period().
    arma::mat draw_period(const VecTvpDiscountPosterior &posterior, arma::uword period,
                          arma::uword draws) const;

    /// Simulates one forecast path per draw from the last in-sample period, in
    /// levels, the coefficients taking one step of their random walk per horizon
    /// under the discount and the cointegration matrix taking none.
    ///
    /// `draws` decides how many, since the posterior it starts from carries
    /// none. `input.forecast.x` is expected in the level layout and not in the
    /// differenced one `input.train.x` uses, as for every VEC here. The level
    /// coefficients are rebuilt from the states at every horizon, because
    /// `A_1 = I + alpha beta' + Gamma_1` is not linear in them.
    ///
    /// The walk is always simulated: `forecast_states` is not read, as for
    /// VarTvpDiscountEstimator, a model whose drift is a declared discount
    /// having said what it does over the horizon already.
    ForecastDraws forecast(const VecTvpDiscountInput &input,
                           const VecTvpDiscountPosterior &posterior, arma::uword draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, 1 x periods rather than draws x periods: the
    /// parameters are integrated out exactly, so there is one number per period
    /// and no Monte Carlo average to take over.
    ///
    /// The density of `/data/train/y`, which for a VEC is in differences. That
    /// is the same number as the density of the level it implies: given the
    /// past, `y_t = y_{t-1} + dy_t` has unit Jacobian, so these are comparable
    /// with a VAR in levels on the same sample as well as with another VEC.
    ///
    /// A consumer that wants the draws x periods matrix WAIC and PSIS-LOO
    /// expect builds it with `draw_period()` instead, and pays for the draws it
    /// then averages over.
    arma::mat log_likelihood(const VecTvpDiscountInput &input,
                             const VecTvpDiscountPosterior &posterior) const;

    /// The log predictive density of what the horizon realised, draws x scored
    /// periods -- one row per draw, one column per period of `input.test.y`,
    /// which for a VEC holds the realised levels.
    ///
    /// **Drawn, where the VAR's is exact, and the reason is the file rather than
    /// the model.** A filter is recursive, so scoring an expanding window costs
    /// nothing in sample and VarTvpDiscountEstimator carries its own recursion
    /// through the realised values. Doing that here would need the error
    /// correction term of each scored period, and a VEC's `/data/forecast/x` is
    /// in the level layout: `w_t` holds the restricted deterministic terms and
    /// the unmodelled variables of the cointegration space, which are not among
    /// the level regressors and cannot be recovered from them. So the score is
    /// simulated, the same way the three sampling VECs score, and the states
    /// step forward exactly as `forecast()` steps them.
    arma::mat predictive_log_density(const VecTvpDiscountInput &input,
                                     const VecTvpDiscountPosterior &posterior,
                                     arma::uword draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_TVP_DISCOUNT_H
