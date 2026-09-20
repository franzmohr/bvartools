// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_discount.h"

#include "core/models/discount_support.h"
#include "core/models/predictive_score.h"

#include <cmath>
#include <stdexcept>

namespace bayests
{

namespace
{

using core::discount_walk_root;
using core::draw_discount_period;
using core::draw_inverse_wishart;
using core::DiscountPath;
using core::DiscountState;
using core::response_by_period;
using core::run_discount_filter;

} // namespace

void VarTvpDiscountInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    if (k == 0)
    {
        throw std::invalid_argument("a model needs at least one endogenous variable");
    }
    if (train.x.n_cols == 0)
    {
        throw std::invalid_argument(
            "the discounted model reads the compact regressors "
            "/data/train/x, not the SUR matrix /data/train/z");
    }

    const arma::uword tt = train.periods(spec.k);
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }
    if (train.x.n_rows != tt)
    {
        throw std::invalid_argument("the regressors and the response disagree on the number of periods");
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

    const arma::uword n_reg = train.x.n_cols;
    if (a_prior.mean.n_rows != n_reg || a_prior.mean.n_cols != k)
    {
        throw std::invalid_argument("the coefficient prior mean must be n_reg x k");
    }
    if (a_prior.cov.n_rows != n_reg || a_prior.cov.n_cols != n_reg)
    {
        throw std::invalid_argument("the coefficient prior covariance must be n_reg x n_reg");
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
    if (spec.delta_sigma < 1.0 &&
        1.0 / (1.0 - spec.delta_sigma) < static_cast<double>(k))
    {
        throw std::invalid_argument(
            "delta_sigma is too small for this many variables: the degrees of "
            "freedom settle at 1 / (1 - delta_sigma), which must exceed k");
    }
}

VarTvpDiscountPosterior
VarTvpDiscountEstimator::estimate(const VarTvpDiscountInput &input,
                                  Reporter &reporter) const
{
    input.validate();

    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword tt = input.train.periods(input.spec.k);

    // The design is the regressors as they stand. A VAR has no block of them
    // that is a function of a parameter, which is the one place a VEC differs
    // and the whole of what VecTvpDiscountEstimator adds.
    DiscountPath path =
        run_discount_filter(response_by_period(input.train, k, tt), input.train.x, input.a_prior,
                            input.u_sigma_prior, input.spec.delta_beta, input.spec.delta_sigma,
                            reporter);

    VarTvpDiscountPosterior posterior;
    posterior.a = std::move(path.a);
    posterior.a_scale = std::move(path.a_scale);
    posterior.a_cov = std::move(path.a_cov);
    posterior.u_sigma = std::move(path.u_sigma);
    posterior.df = std::move(path.df);
    posterior.forecast_mean = std::move(path.forecast_mean);
    posterior.loglik = std::move(path.loglik);

    reporter.finish();
    return posterior;
}

arma::mat VarTvpDiscountEstimator::draw_period(const VarTvpDiscountPosterior &posterior,
                                               const arma::uword period,
                                               const arma::uword draws) const
{
    const arma::uword k = static_cast<arma::uword>(
        std::lround(std::sqrt(static_cast<double>(posterior.u_sigma.n_rows))));

    return draw_discount_period(posterior.a, posterior.a_cov, posterior.u_sigma, posterior.df,
                                period, draws, k);
}

arma::mat VarTvpDiscountEstimator::log_likelihood(const VarTvpDiscountInput &input,
                                                  const VarTvpDiscountPosterior &posterior) const
{
    // One row, not one per draw: the parameters are integrated out exactly, so
    // there is nothing to average over. The caller that wants the draws x
    // periods matrix WAIC and PSIS-LOO expect builds it from draw_period().
    (void)input;
    return arma::mat(posterior.loglik.t());
}

ForecastDraws VarTvpDiscountEstimator::forecast(const VarTvpDiscountInput &input,
                                                const VarTvpDiscountPosterior &posterior,
                                                const arma::uword draws,
                                                Reporter &reporter) const
{
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword h = static_cast<arma::uword>(input.spec.h);
    if (h == 0)
    {
        throw std::invalid_argument("no forecast horizon was requested");
    }
    if (input.forecast.x.n_rows != h)
    {
        throw std::invalid_argument("the forecast regressors do not cover the horizon");
    }

    const arma::uword last = posterior.periods() - 1;
    const arma::uword n_reg = input.forecast.x.n_cols;

    const arma::mat c = arma::reshape(posterior.a_cov.col(last), n_reg, n_reg);
    const arma::mat s = arma::reshape(posterior.u_sigma.col(last), k, k);
    const double df = posterior.df(last);

    // The random walk innovation the discount implies, held at the last period
    // the sample pinned down.
    const arma::mat chol_step = discount_walk_root(c, input.spec.delta_beta);

    const arma::mat start = draw_period(posterior, last, draws);

    arma::mat chol_inv_scale;
    if (!arma::chol(chol_inv_scale, arma::inv_sympd(arma::symmatu(df * s)), "lower"))
    {
        throw std::invalid_argument("the final posterior scale is not positive definite");
    }

    ForecastDraws out;
    out.values.set_size(h * k, draws);

    for (arma::uword i = 0; i < draws; i++)
    {
        reporter.check_interrupt();

        arma::mat theta = arma::reshape(start.col(i), k, n_reg).t();
        const arma::mat sigma = draw_inverse_wishart(chol_inv_scale, df, k);
        arma::mat chol_sigma;
        if (!arma::chol(chol_sigma, arma::symmatu(sigma), "lower"))
        {
            throw std::invalid_argument("a drawn error covariance was not positive definite");
        }

        arma::mat regressors = input.forecast.x;
        for (arma::uword step = 0; step < h; step++)
        {
            theta += chol_step * arma::mat(n_reg, k, arma::fill::randn);

            const arma::vec z = regressors.row(step).t();
            const arma::vec value = theta.t() * z + chol_sigma * arma::vec(k, arma::fill::randn);
            out.values.submat(step * k, i, (step + 1) * k - 1, i) = value;

            // Carry the realisation into the lag blocks of the next period, the
            // way every other model in this project simulates a path.
            if (step + 1 < h)
            {
                const arma::uword lags = static_cast<arma::uword>(input.spec.p);
                for (arma::uword j = lags; j-- > 1;)
                {
                    regressors.submat(step + 1, j * k, step + 1, (j + 1) * k - 1) =
                        regressors.submat(step, (j - 1) * k, step, j * k - 1);
                }
                if (lags > 0)
                {
                    regressors.submat(step + 1, 0, step + 1, k - 1) = value.t();
                }
            }
        }

        reporter.progress(static_cast<long long>(i) + 1,
                          static_cast<long long>(draws));
    }

    reporter.finish();
    return out;
}

arma::mat
VarTvpDiscountEstimator::predictive_log_density(const VarTvpDiscountInput &input,
                                                const VarTvpDiscountPosterior &posterior) const
{
    if (input.test.y.n_elem == 0)
    {
        throw std::invalid_argument("the model file carries no realised values to score against");
    }

    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword scored = input.test.y.n_rows;

    // The lag blocks of `/data/forecast/x` hold whatever the caller put there,
    // which for a host that has only asked for a forecast is a placeholder: a
    // forecast overwrites them as it simulates. This recursion does not
    // simulate, so it has to fill them from what was realised, exactly as every
    // other VAR here scores. Reading the rows raw scored the first period
    // correctly and then fed the filter a placeholder, which left the state --
    // and every period after the first -- not a number.
    const arma::mat x =
        core::realised_regressors(input.forecast.x, input.test.y, input.spec.k, input.spec.p);

    const arma::uword last = posterior.periods() - 1;
    const arma::uword n_reg = input.forecast.x.n_cols;

    // Carrying the filter through the realised values, which is the same
    // recursion estimate() ran and needs no re-estimation and no draws: each
    // period conditions on everything realised before it by construction.
    DiscountState state;
    state.m = arma::reshape(posterior.a.col(last), k, n_reg).t();
    state.c = arma::reshape(posterior.a_cov.col(last), n_reg, n_reg);
    state.df = posterior.df(last);
    state.d = state.df * arma::reshape(posterior.u_sigma.col(last), k, k);

    arma::mat out(1, scored);
    for (arma::uword t = 0; t < scored; t++)
    {
        out(0, t) = state.step(x.row(t).t(), input.test.y.row(t).t(),
                               input.spec.delta_beta, input.spec.delta_sigma);
    }
    return out;
}

} // namespace bayests
