// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_tvp_discount.h"

#include "core/models/discount_support.h"
#include "core/models/vec_support.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace bayests
{

namespace
{

using core::discount_walk_root;
using core::draw_discount_period;
using core::draw_inverse_wishart;
using core::DiscountPath;
using core::response_by_period;
using core::run_discount_filter;
using core::score_vec_forecast;
using core::simulate_vec_forecast;
using core::VecForecastStep;

/// The design the filter runs against: the `rank` error correction columns in
/// front of the compact regressors, one row per period.
///
/// `beta' w_t` is the t-th row of `w beta`, so the whole block is one product
/// and the loop every VEC sampler runs to rebuild it per draw is not needed --
/// the space does not move, and this is built once for the whole sample. That
/// single line is the entire difference between this model and VarTvpDiscount.
arma::mat build_design(const VecTvpDiscountInput &input, const arma::uword tt)
{
    const arma::uword rank = static_cast<arma::uword>(input.spec.rank);
    const arma::uword n_x = static_cast<arma::uword>(input.spec.n_x_vec());

    arma::mat design(tt, rank + n_x);
    if (rank > 0)
    {
        design.cols(0, rank - 1) = input.train.w * input.beta;
    }
    if (n_x > 0)
    {
        design.cols(rank, rank + n_x - 1) = input.train.x;
    }
    return design;
}

/// How a VecTvpDiscount's states move over a horizon, in one place.
///
/// The forecast and the score need the same walk and have to take it
/// identically, or the two would describe different models from the same file:
/// the same reason VecTvpWishartWalk exists, and called the same way, once per
/// draw and horizon with i = 0 first.
///
/// What moves is the coefficient matrix, by the random walk `delta_beta`
/// implies at the last in-sample period, and nothing else. The cointegration
/// matrix is handed over unchanged at every horizon, which is the model; the
/// error covariance is drawn once per draw and held, the discounted Wishart
/// having no innovation to step and no mean to revert to.
struct VecTvpDiscountWalk
{
    const VecTvpDiscountPosterior &posterior;
    arma::uword k;
    arma::uword n_design;
    bool use_beta;

    arma::mat start;          ///< nparams x draws, the last period drawn from.
    arma::mat chol_step;      ///< n_design square, the walk's innovation factor.
    arma::mat chol_inv_scale; ///< k square, what Sigma is drawn through.
    double df = 0.0;

    arma::mat theta; ///< n_design x k, where the draw being simulated has got to.

    VecTvpDiscountWalk(const VecTvpDiscountInput &input,
                       const VecTvpDiscountPosterior &posterior_in, const arma::uword draws,
                       const VecTvpDiscountEstimator &estimator)
        : posterior(posterior_in), k(static_cast<arma::uword>(input.spec.k)),
          n_design(static_cast<arma::uword>(input.n_design())), use_beta(input.use_beta())
    {
        if (posterior.periods() == 0)
        {
            throw std::invalid_argument("the posterior carries no periods to forecast from");
        }
        const arma::uword last = posterior.periods() - 1;
        if (posterior.a.n_rows != n_design * k)
        {
            throw std::invalid_argument(
                "the posterior and the model disagree on the width of a period: `a` has " +
                std::to_string(posterior.a.n_rows) + " rows against the " +
                std::to_string(n_design * k) + " this model carries");
        }
        if (use_beta && posterior.beta.n_elem != static_cast<arma::uword>(input.spec.n_beta()))
        {
            throw std::invalid_argument(
                "the posterior does not carry the cointegration matrix the forecast is taken "
                "over; a discounted VEC stores the space it conditioned on beside its draws");
        }

        df = posterior.df(last);
        start = estimator.draw_period(posterior, last, draws);
        chol_step = discount_walk_root(
            arma::reshape(posterior.a_cov.col(last), n_design, n_design), input.spec.delta_beta);

        const arma::mat s = arma::reshape(posterior.u_sigma.col(last), k, k);
        if (!arma::chol(chol_inv_scale, arma::inv_sympd(arma::symmatu(df * s)), "lower"))
        {
            throw std::invalid_argument("the final posterior scale is not positive definite");
        }
    }

    void operator()(const arma::uword draw, const int i, VecForecastStep &out)
    {
        if (i == 0)
        {
            theta = arma::reshape(start.col(draw), k, n_design).t();

            const arma::mat sigma = draw_inverse_wishart(chol_inv_scale, df, k);
            out.period.u_sigma_inv = arma::vectorise(arma::inv_sympd(arma::symmatu(sigma)));
            if (!arma::chol(out.error_root, arma::symmatu(sigma), "lower"))
            {
                throw std::invalid_argument("a drawn error covariance was not positive definite");
            }
            if (use_beta)
            {
                out.period.beta = posterior.beta;
            }
        }

        // The state at step i is the one the observation at step i is generated
        // from, so the first innovation lands before the first horizon. See
        // core/models/forecast_states.h.
        theta += chol_step * arma::mat(n_design, k, arma::fill::randn);
        out.period.a = arma::vectorise(theta.t());
    }
};

} // namespace

void VecTvpDiscountInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    if (k == 0)
    {
        throw std::invalid_argument("a model needs at least one endogenous variable");
    }
    if (spec.rank < 0 || spec.k_beta < 0)
    {
        throw std::invalid_argument("the cointegration rank and k_beta cannot be negative");
    }
    if (n_design() == 0)
    {
        throw std::invalid_argument(
            "the discounted VEC has no regressors at all: with a rank of zero it needs the "
            "compact regressors /data/train/x, which it reads in place of the SUR matrix "
            "/data/train/z");
    }

    const arma::uword tt = train.periods(spec.k);
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    const arma::uword rank = static_cast<arma::uword>(spec.rank);
    const arma::uword k_beta = static_cast<arma::uword>(spec.k_beta);
    const arma::uword n_x = static_cast<arma::uword>(spec.n_x_vec());

    if (n_x > 0 && (train.x.n_rows != tt || train.x.n_cols != n_x))
    {
        throw std::invalid_argument(
            "the compact regressors must be tt x " + std::to_string(n_x) + ", got " +
            std::to_string(train.x.n_rows) + " x " + std::to_string(train.x.n_cols) +
            "; the discounted VEC reads /data/train/x and not the SUR matrix /data/train/z");
    }

    if (use_beta())
    {
        if (train.w.n_rows != tt || train.w.n_cols != k_beta)
        {
            throw std::invalid_argument(
                "the error correction term must be tt x k_beta, that is " + std::to_string(tt) +
                " x " + std::to_string(k_beta) + ", got " + std::to_string(train.w.n_rows) + " x " +
                std::to_string(train.w.n_cols));
        }
        // The one input this model has and its sampling siblings do not. They
        // draw the space; this conditions on it, so a missing or misshapen beta
        // is a missing model rather than a missing starting value.
        if (beta.n_rows != k_beta || beta.n_cols != rank)
        {
            throw std::invalid_argument(
                "the discounted VEC conditions on a fixed cointegration matrix, which must be "
                "k_beta x rank, that is " + std::to_string(k_beta) + " x " + std::to_string(rank) +
                ", got " + std::to_string(beta.n_rows) + " x " + std::to_string(beta.n_cols));
        }
        if (!beta.is_finite())
        {
            throw std::invalid_argument("the cointegration matrix holds a value that is not finite");
        }
    }

    // Both discounts are models in their own right at one, so the range is
    // half open rather than open at that end.
    if (!(spec.delta_beta > 0.0 && spec.delta_beta <= 1.0))
    {
        throw std::invalid_argument("delta_beta must lie in (0, 1]");
    }
    if (!(spec.delta_sigma > 0.0 && spec.delta_sigma <= 1.0))
    {
        throw std::invalid_argument("delta_sigma must lie in (0, 1]");
    }

    if (spec.varsel != VarSelection::none)
    {
        throw std::invalid_argument(
            "variable selection is not available for the discounted model: it "
            "has no draws for an inclusion indicator to be drawn alongside");
    }

    // Nothing is iterated, so there is nothing to discard and nothing to keep
    // one of. Refused rather than ignored because both reach the file: the
    // `start`, `end` and `thin` attributes written beside a forecast come from
    // `thin`, and a file that asked for a burn-in and got draws labelled as
    // though it had happened is output that looks like output. `iterations`
    // keeps its meaning -- how many i.i.d. draws a stage that needs them takes.
    if (spec.burnin != 0)
    {
        throw std::invalid_argument(
            "the discounted model has no chain to burn in: burnin must be zero, and "
            "iterations alone says how many i.i.d. draws a forecast takes");
    }
    if (spec.thin != 1)
    {
        throw std::invalid_argument(
            "the discounted model draws i.i.d. rather than sweeping, so there is nothing to "
            "thin: thin must be one");
    }
    if (spec.structural)
    {
        throw std::invalid_argument(
            "a structural model cannot be estimated against an unrestricted "
            "error covariance, which an inverse Wishart posterior is");
    }

    const arma::uword n_design_u = static_cast<arma::uword>(n_design());
    if (a_prior.mean.n_rows != n_design_u || a_prior.mean.n_cols != k)
    {
        throw std::invalid_argument(
            "the coefficient prior mean must be n_design x k, that is " +
            std::to_string(n_design_u) + " x " + std::to_string(k) + ", got " +
            std::to_string(a_prior.mean.n_rows) + " x " + std::to_string(a_prior.mean.n_cols));
    }
    if (a_prior.cov.n_rows != n_design_u || a_prior.cov.n_cols != n_design_u)
    {
        throw std::invalid_argument("the coefficient prior covariance must be n_design x n_design");
    }
    if (u_sigma_prior.df <= 0)
    {
        throw std::invalid_argument("the prior degrees of freedom must be positive");
    }
    if (u_sigma_prior.scale.n_rows != k || u_sigma_prior.scale.n_cols != k)
    {
        throw std::invalid_argument("the error covariance prior scale must be k x k");
    }

    // The degrees of freedom settle at 1 / (1 - delta_sigma) whatever the prior
    // says, so a short memory cannot carry a wide system. Refused rather than
    // warned about: below k the inverse Wishart is improper and every scale
    // read off it is meaningless.
    if (spec.delta_sigma < 1.0 && 1.0 / (1.0 - spec.delta_sigma) < static_cast<double>(k))
    {
        throw std::invalid_argument(
            "delta_sigma is too small for this many variables: the degrees of "
            "freedom settle at 1 / (1 - delta_sigma), which must exceed k");
    }
}

VecTvpDiscountPosterior VecTvpDiscountEstimator::estimate(const VecTvpDiscountInput &input,
                                                          Reporter &reporter) const
{
    input.validate();

    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword tt = input.train.periods(input.spec.k);

    DiscountPath path = run_discount_filter(response_by_period(input.train, k, tt),
                                            build_design(input, tt), input.a_prior,
                                            input.u_sigma_prior, input.spec.delta_beta,
                                            input.spec.delta_sigma, reporter);

    VecTvpDiscountPosterior posterior;
    posterior.a = std::move(path.a);
    posterior.a_scale = std::move(path.a_scale);
    posterior.a_cov = std::move(path.a_cov);
    posterior.u_sigma = std::move(path.u_sigma);
    posterior.df = std::move(path.df);
    posterior.forecast_mean = std::move(path.forecast_mean);
    posterior.loglik = std::move(path.loglik);

    // Carried through rather than estimated, so that what comes out says which
    // space it conditioned on. See VecTvpDiscountPosterior::beta.
    if (input.use_beta())
    {
        posterior.beta = arma::vectorise(input.beta);
    }

    reporter.finish();
    return posterior;
}

arma::mat VecTvpDiscountEstimator::draw_period(const VecTvpDiscountPosterior &posterior,
                                               const arma::uword period,
                                               const arma::uword draws) const
{
    const arma::uword k = static_cast<arma::uword>(
        std::lround(std::sqrt(static_cast<double>(posterior.u_sigma.n_rows))));

    return draw_discount_period(posterior.a, posterior.a_cov, posterior.u_sigma, posterior.df,
                                period, draws, k);
}

arma::mat VecTvpDiscountEstimator::log_likelihood(const VecTvpDiscountInput &input,
                                                  const VecTvpDiscountPosterior &posterior) const
{
    // One row, not one per draw: the parameters are integrated out exactly, so
    // there is nothing to average over. See the header for why this is
    // comparable with a model estimated in levels.
    (void)input;
    return arma::mat(posterior.loglik.t());
}

ForecastDraws VecTvpDiscountEstimator::forecast(const VecTvpDiscountInput &input,
                                                const VecTvpDiscountPosterior &posterior,
                                                const arma::uword draws,
                                                Reporter &reporter) const
{
    if (draws == 0)
    {
        throw std::invalid_argument("a draw count of zero forecasts nothing");
    }

    VecTvpDiscountWalk walk(input, posterior, draws, *this);
    return ForecastDraws{
        simulate_vec_forecast(input.spec, input.forecast, draws, reporter, walk)};
}

arma::mat VecTvpDiscountEstimator::predictive_log_density(const VecTvpDiscountInput &input,
                                                          const VecTvpDiscountPosterior &posterior,
                                                          const arma::uword draws) const
{
    if (input.test.y.n_elem == 0)
    {
        throw std::invalid_argument("the model file carries no realised values to score against");
    }
    if (draws == 0)
    {
        throw std::invalid_argument("a draw count of zero scores nothing");
    }

    VecTvpDiscountWalk walk(input, posterior, draws, *this);
    return score_vec_forecast(input.spec, input.forecast, input.test.y, draws, walk);
}

} // namespace bayests
