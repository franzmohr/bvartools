// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_wishart.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/wishart.h"
#include "core/models/model_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_normal_precision;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarTvpWishartDraws VarTvpWishartSampler::draw_coefficients(const VarTvpWishartInput &input,
                                                       Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int kk = k * k;
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;

    const int nparams = static_cast<int>(z.n_cols);
    const bool use_a = nparams > 0;
    const int tt = static_cast<int>(y.n_elem) / k;
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    const arma::mat diag_tt = arma::eye<arma::mat>(tt, tt);
    arma::mat ymat = arma::reshape(y, k, tt);

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VarTvpWishartDraws out;

    // Coefficients: a whole path per draw, moved by the simulation smoother.
    arma::mat a, a_B, a_sigma, a_lag;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    const arma::vec &a_sigma_prior_rate = input.a_prior.sigma.rate;

    arma::vec a0;
    arma::mat a0_post_v, a0_sigma_inv;
    const arma::vec &a0_prior_mu = input.a_prior.initial_state.mu;
    const arma::mat &a0_prior_v_inv = input.a_prior.initial_state.v_inv;

    // Variable selection
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;
    arma::vec a_theta_res;

    if (use_a)
    {
        a = input.initial.a;
        a_lag = a;
        a_sigma = input.initial.a_sigma_inv;
        a_sigma.diag() = 1 / a_sigma.diag();
        out.a = arma::mat(nparams * tt, iterations);
        out.a_sigma = arma::mat(nparams, iterations);
        a_B = arma::eye<arma::mat>(nparams, nparams);

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
                a_theta_res = arma::zeros<arma::vec>(y.n_elem);
            }
        }
    }


    // Error term
    arma::mat sse;
    const int post_u_sigma_df = input.u_sigma_prior.df + tt;
    const arma::mat &prior_u_sigma_scale = input.u_sigma_prior.scale;
    arma::mat u_sigma_inv = input.initial.u_sigma_inv;
    arma::mat u = arma::reshape(y, k, tt);
    out.u_sigma_inv = arma::mat(kk, iterations);

    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);

    arma::sp_mat u_sigma_inv_diag;
    const arma::sp_mat diag_tt_sp = arma::speye<arma::sp_mat>(tt, tt);


    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            // The precision is the same in every period, so the block diagonal
            // is kron(I_tt, u_sigma_inv) and the smoother's covariance is one
            // k x k inverse repeated. Both are rebuilt here rather than before
            // the loop because u_sigma_inv is redrawn at the end of every
            // iteration: built once, the smoother would condition on the
            // initial covariance for the whole chain.
            //
            // The inverse is a dense solve of a k x k matrix. spsolve() would
            // need SuperLU, which an embedded host such as an R package does
            // not have, and there is nothing sparse about a block this size.
            u_sigma = arma::repmat(arma::solve(u_sigma_inv, diag_k), tt, 1);
            u_sigma_inv_diag = arma::kron(diag_tt_sp, arma::sp_mat(u_sigma_inv));

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a
            a = kalman_durbin_koopman_2002(ymat, z, u_sigma, a_sigma, a_B, a0, a_sigma)
                    .cols(0, tt - 1);

            // Draw a_sigma
            a_lag.col(0) = a0;
            a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
            a_lag = a - a_lag;
            a_sigma_post_scale = 1 / (a_sigma_prior_rate + arma::sum(arma::pow(a_lag, 2), 1) * 0.5);
            for (int i = 0; i < nparams; i++)
            {
                a_sigma(i, i) = 1 / arma::randg<double>(arma::distr_param(a_sigma_post_shape(i), a_sigma_post_scale(i)));
            }

            // Draw a0
            a0_sigma_inv.diag() = 1 / a_sigma.diag();
            a0_post_v = a0_prior_v_inv + a0_sigma_inv;
            a0 = draw_normal_precision(a0_post_v,
                                       a0_prior_v_inv * a0_prior_mu + a0_sigma_inv * a.col(0));

            if (a_bvs)
            {
                z = z_bvs;
                bvs_sweep(*a_bvs, a, BvsScope::path_row, [&](const arma::mat &theta) {
                    for (int i = 0; i < tt; i++)
                    {
                        a_theta_res.subvec(i * k, (i + 1) * k - 1) =
                            y.subvec(i * k, (i + 1) * k - 1) -
                            z.rows(i * k, (i + 1) * k - 1) * theta.col(i);
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

        // Update u_sigma_inv
        u_sigma_inv = wishart(u, prior_u_sigma_scale, post_u_sigma_df);
        

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            // a
            if (use_a)
            {
                out.a.col(draw_pos) = arma::vectorise(a);
                out.a_sigma.col(draw_pos) = arma::vectorise(a_sigma.diag());
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_bvs->lambda;
                }
            }

            // Measurement error
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarTvpWishartSampler::forecast(const VarTvpWishartInput &input,
                                           const VarTvpWishartDraws &coefficients,
                                           Reporter &reporter) const
{
    const int k = input.spec.k;
    const int p = input.spec.p;
    const int h = input.spec.h;
    const bool structural = input.spec.structural;
    const int n_structural = input.spec.n_structural();
    const int n_non_structural = input.spec.n_non_structural();
    const int nparams = n_non_structural + n_structural;

    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (h <= 0)
    {
        throw std::invalid_argument("forecast horizon (h) must be positive");
    }
    if (coefficients.u_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
    }

    arma::mat z = input.forecast.z;

    require_forecast_regressors(input.spec, z);

    // Counted off the model's dimensions rather than off `z`: the coefficients
    // move with time, so what the forecast starts from is the last in-sample
    // period of the posterior, and the contemporaneous block at the end of it
    // has no column in `z` to be counted by.
    const bool use_a = n_non_structural > 0;

    if (nparams > 0 && !coefficients.has_a())
    {
        throw std::invalid_argument("forecasting this model needs posterior draws of a, which "
                                    "are missing");
    }

    // The caller hands over the period the forecast starts from, as the header
    // says; slicing the path here as well would take the last period of a
    // matrix that is already one period wide.
    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

    if (use_a && z.n_cols != a.n_rows)
    {
        throw std::invalid_argument(
            "forecast regressors and coefficient draws disagree: z has " +
            std::to_string(z.n_cols) + " columns, a has " + std::to_string(a.n_rows) +
            " rows after the structural split");
    }

    const arma::uword draws = coefficients.iterations();
    const bool p_larger_than_0 = p > 0;

    arma::mat fcst = arma::zeros<arma::mat>(h * k, draws);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::vec eigval;
    arma::mat eigvec;

    // Calculate forecasts
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        // Once per draw: nothing in it depends on the horizon.
        const arma::mat a0_inv =
            structural ? structural_inverse(a0, draw, diag_k) : arma::mat();

        for (int i = 0; i < h; i++)
        {
            if (use_a)
            {
                // Update z
                if (i > 0 && p_larger_than_0)
                {
                    update_forecast_lags(z, fcst, draw, i, k, p, diag_k);
                }
                // Update forecast
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = z.rows(i * k, (i + 1) * k - 1) * a.col(draw);
            }

            // Add error
            arma::eig_sym(eigval, eigvec, arma::solve(arma::reshape(coefficients.u_sigma_inv.col(draw), k, k), diag_k));
            fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::trans(eigvec) * arma::randn(k);

            // A_0 y_t = A_1 y_{t-1} + ... + u_t, so the inverse applies to the
            // whole right-hand side, signal and error alike.
            if (structural)
            {
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) =
                    a0_inv * fcst.submat(i * k, draw, (i + 1) * k - 1, draw);
            }
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat VarTvpWishartSampler::log_likelihood(const VarTvpWishartInput &input,
                                               const VarTvpWishartDraws &coefficients) const
{
    const int k = input.spec.k;

    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (coefficients.u_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
    }

    const arma::vec y = stacked_response(input.train);
    const arma::mat &z = input.train.z;
    const int nparams = static_cast<int>(z.n_cols);
    const bool use_a = nparams > 0;

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }

    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    arma::mat loglik = arma::mat(draws, tt);

    // Calculate errors. Each period has its own coefficients, so the design is
    // the block diagonal of the per-period regressors rather than `z` itself.
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        arma::sp_mat z_large = arma::zeros<arma::sp_mat>(k * tt, nparams * tt);
        for (int i = 0; i < tt; i++)
        {
            z_large.submat(i * k, i * nparams, (i + 1) * k - 1, (i + 1) * nparams - 1) = z.rows(i * k, (i + 1) * k - 1);
        }
        u = u - z_large * coefficients.a;
    }

    // Calculate log likelihood
    const arma::mat diag_k = arma::eye(k, k);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, k);
        const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
        for (int i = 0; i < tt; i++)
        {
            const double part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
