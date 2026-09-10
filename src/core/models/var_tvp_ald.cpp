// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_ald.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/models/ald_support.h"
#include "core/models/model_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::AldShape;
using core::ald_log_density;
using core::ald_shape;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_ald_scale;
using core::draw_ald_weights;
using core::draw_normal_precision;
using core::stacked_response;

VarTvpAldDraws VarTvpAldSampler::draw_coefficients(const VarTvpAldInput &input,
                                                   Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;

    const int nparams = static_cast<int>(z.n_cols);
    const bool use_a = nparams > 0;
    const int tt = static_cast<int>(y.n_elem) / k;

    const arma::mat ymat = arma::reshape(y, k, tt);

    const AldShape shape = ald_shape(input.spec.quantile);

    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VarTvpAldDraws out;

    // Coefficients
    arma::mat a, a_B, a_sigma, a_lag;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    arma::vec a0;
    arma::mat a0_post_v, a0_sigma_inv;

    // Variable selection
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;
    arma::vec a_theta_res;

    const arma::vec &a_sigma_prior_rate = input.a_prior.sigma.rate;
    const arma::vec &a0_prior_mu = input.a_prior.initial_state.mu;
    const arma::mat &a0_prior_v_inv = input.a_prior.initial_state.v_inv;

    if (use_a)
    {
        a = input.initial.a;
        a_lag = a;
        a_sigma = input.initial.a_sigma_inv;
        a_sigma.diag() = 1 / a_sigma.diag();
        a_B = arma::eye<arma::mat>(nparams, nparams);

        out.a = arma::mat(nparams * tt, iterations);
        out.a_sigma = arma::mat(nparams, iterations);

        a_sigma_post_shape = input.a_prior.sigma.shape + tt * 0.5;
        a_sigma_post_scale = a_sigma_prior_rate;

        a0 = input.initial.a_init;
        a0_sigma_inv = a_sigma;
        a0_sigma_inv.diag() = 1 / a_sigma.diag();

        if (use_varsel)
        {
            out.a_lambda = arma::mat(nparams, iterations);

            if (use_bvs)
            {
                z_bvs = z;
                a_bvs.emplace(input.initial.a_lambda, input.a_varsel_prior);
                a_theta_res = arma::zeros<arma::vec>(k * tt);
            }
        }
    }

    // The error term and the latent scales it is a mixture over.
    arma::mat u = ymat;
    arma::mat w = input.initial.w;
    arma::vec u_scale = input.initial.u_scale;

    // What the smoother measures against: the data less the skew the mixture
    // carries, k x tt to match `ymat`.
    arma::mat y_adjusted = ymat - shape.theta * arma::trans(w);

    // The measurement covariance the smoother wants, as a (k tt) x k stack of
    // per-period blocks. Every block is diagonal -- this model has no covariance
    // block -- so it is filled rather than inverted, which is the one place this
    // is cheaper than the stochastic volatility model it is otherwise a copy of.
    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);

    // The same information as a precision, which is what the BVS sweep scores
    // against.
    arma::sp_mat u_sigma_inv_diag = arma::eye<arma::sp_mat>(k * tt, k * tt);

    const arma::vec u_scale_post_shape =
        input.u_scale_prior.shape + static_cast<double>(tt) * 3.0 / 2.0;
    const arma::vec &u_scale_prior_rate = input.u_scale_prior.rate;

    // Fills both representations from the current latent scales.
    const auto rebuild_variances = [&]() {
        for (int t = 0; t < tt; t++)
        {
            for (int i = 0; i < k; i++)
            {
                const double variance = shape.tau2 * u_scale(i) * w(t, i);
                u_sigma(k * t + i, i) = variance;
                u_sigma_inv_diag(k * t + i, k * t + i) = 1.0 / variance;
            }
        }
    };
    rebuild_variances();

    out.u_scale = arma::mat(k, iterations);
    out.u_omega_inv = arma::mat(k * tt, iterations);
    out.u_sigma_inv = arma::mat(k * k * tt, iterations);

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a. The response carries the offset: the smoother measures
            // z_t a_t against y_t - theta w_t, not against y_t.
            a = kalman_durbin_koopman_2002(y_adjusted, z, u_sigma, a_sigma, a_B, a0, a_sigma)
                    .cols(0, tt - 1);

            // Draw a_sigma
            a_lag.col(0) = a0;
            a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
            a_lag = a - a_lag;
            a_sigma_post_scale = 1 / (a_sigma_prior_rate + arma::sum(arma::pow(a_lag, 2), 1) * 0.5);
            for (int i = 0; i < nparams; i++)
            {
                a_sigma(i, i) =
                    1 / arma::randg<double>(
                            arma::distr_param(a_sigma_post_shape(i), a_sigma_post_scale(i)));
            }

            // Draw a0
            a0_sigma_inv.diag() = 1 / a_sigma.diag();
            a0_post_v = a0_prior_v_inv + a0_sigma_inv;
            a0 = draw_normal_precision(a0_post_v,
                                       a0_prior_v_inv * a0_prior_mu + a0_sigma_inv * a.col(0));

            if (a_bvs)
            {
                z = z_bvs;
                // path_row rather than element, for the reason spelled out at
                // the same call in var_tvp_gamma.cpp: element scope reached
                // period 0 alone. The residual carries the offset here too.
                bvs_sweep(*a_bvs, a, BvsScope::path_row, [&](const arma::mat &theta) {
                    for (int i = 0; i < tt; i++)
                    {
                        a_theta_res.subvec(i * k, (i + 1) * k - 1) =
                            y_adjusted.col(i) - z.rows(i * k, (i + 1) * k - 1) * theta.col(i);
                    }
                    return -arma::as_scalar(arma::trans(a_theta_res) * u_sigma_inv_diag *
                                            a_theta_res) /
                           2;
                });
            }

            for (int i = 0; i < tt; i++)
            {
                u.col(i) = ymat.col(i) - z.rows(i * k, (i + 1) * k - 1) * a.col(i);
            }
        }
        else
        {
            u = ymat;
        }

        // Update the latent scales, one per observation.
        w = draw_ald_weights(u, u_scale, shape);

        // Update the scale of the asymmetric Laplace, one per equation.
        u_scale = draw_ald_scale(u, w, shape, u_scale_post_shape, u_scale_prior_rate);

        // Rebuild what the next sweep measures and weights against.
        y_adjusted = ymat - shape.theta * arma::trans(w);
        rebuild_variances();

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            if (use_a)
            {
                out.a.col(draw_pos) = arma::vectorise(a);
                out.a_sigma.col(draw_pos) = a_sigma.diag();
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_bvs->lambda;
                }
            }

            out.u_scale.col(draw_pos) = u_scale;
            out.u_omega_inv.col(draw_pos) = u_sigma_inv_diag.diag();

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * k * k, draw_pos, (i + 1) * k * k - 1, draw_pos) =
                    arma::vectorise(arma::mat(u_sigma_inv_diag.submat(
                        i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1)));
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarTvpAldSampler::forecast(const VarTvpAldInput &input, const VarTvpAldDraws &draws,
                                         Reporter &reporter) const
{
    (void)input;
    (void)draws;
    (void)reporter;

    throw std::invalid_argument(
        "a quantile regression model does not forecast: the h step ahead quantile is not the "
        "quantile of the iterated one step ahead quantiles, so iterating this model forward "
        "would produce a path that cannot be read as a quantile of anything");
}

arma::mat VarTvpAldSampler::log_likelihood(const VarTvpAldInput &input,
                                           const VarTvpAldDraws &coefficients) const
{
    const int k = input.spec.k;

    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (coefficients.u_scale.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of the asymmetric Laplace scale are missing");
    }

    const arma::vec y = stacked_response(input.train);
    const arma::mat &z = input.train.z;
    const bool use_a = z.n_cols > 0;

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }

    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    const int nparams = static_cast<int>(z.n_cols);
    const double q = input.spec.quantile;

    arma::mat loglik(draws, tt);
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::vec scale = coefficients.u_scale.col(draw);

        // One coefficient vector per period, so the residual is formed period by
        // period rather than in one product.
        for (int i = 0; i < tt; i++)
        {
            arma::vec residual = y.subvec(i * k, (i + 1) * k - 1);
            if (use_a)
            {
                const arma::vec a_i =
                    coefficients.a.submat(i * nparams, draw, (i + 1) * nparams - 1, draw);
                residual -= z.rows(i * k, (i + 1) * k - 1) * a_i;
            }

            double total = 0.0;
            for (int j = 0; j < k; j++)
            {
                total += ald_log_density(residual(j), scale(j), q);
            }
            loglik(draw, static_cast<arma::uword>(i)) = total;
        }
    }

    return loglik;
}

} // namespace bayests
