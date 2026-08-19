// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_stochvol.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/model_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::build_psi_regressors;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_normal_precision;
using core::fill_psi_path;
using core::fill_strict_lower_triangle;
using core::stacked_response;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarTvpStochvolDraws VarTvpStochvolSampler::draw_coefficients(const VarTvpStochvolInput &input,
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
    const arma::mat ymat = arma::reshape(y, k, tt);

    const bool use_psi = input.use_psi();
    const int n_psi = input.spec.n_psi();

    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;
    const bool use_varsel_psi = use_psi && input.psi_varsel == VarSelection::bvs;

    VarTvpStochvolDraws out;

    // Coefficients
    arma::mat a, a_B, a_sigma, a_lag, post_a_v;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    arma::vec a0, a0_post_mu;
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

    // Covariance block
    arma::mat psi, Psi, psi_B, psi_lag, psi_sigma, psi_u_omega, psi_y, psi_z, psi_z_bvs;
    arma::vec psi_sigma_post_shape, psi_sigma_post_scale;
    arma::vec psi0, psi0_post_mu;
    arma::mat psi0_post_v, psi0_sigma_inv;

    std::optional<BvsBlock> psi_bvs;
    arma::vec psi_theta_res;
    arma::mat Psi_lambda, psi_u_omega_inv_diag;

    if (use_psi)
    {
        psi = input.initial.psi;
        psi_lag = psi;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        psi_B = arma::eye<arma::mat>(n_psi, n_psi);
        Psi = arma::eye<arma::mat>(k * tt, k * tt);
        psi_u_omega = arma::zeros<arma::mat>((k - 1) * tt, k - 1);

        out.psi = arma::mat(kk * tt, iterations);
        out.psi_sigma = arma::mat(n_psi, iterations);

        fill_psi_path(Psi, psi, k);

        psi_sigma = input.initial.psi_sigma_inv;
        psi_sigma.diag() = 1 / psi_sigma.diag();
        psi_sigma_post_shape = input.psi_prior.sigma.shape + tt * 0.5;
        psi_sigma_post_scale = input.psi_prior.sigma.rate;

        psi0 = input.initial.psi_init;
        psi0_sigma_inv = psi_sigma;
        psi0_sigma_inv.diag() = 1 / psi_sigma.diag();

        if (use_varsel_psi)
        {
            out.psi_lambda = arma::mat(kk, iterations);
            Psi_lambda = arma::eye<arma::mat>(k, k);
            psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            psi_theta_res = arma::zeros<arma::vec>((k - 1) * tt);
            psi_u_omega_inv_diag = arma::zeros<arma::mat>((k - 1) * tt, (k - 1) * tt);
        }
    }

    // Error variances
    arma::mat u = arma::reshape(y, k, tt);
    const arma::vec &h_y_offset = input.u_sigma_prior.offset;

    arma::vec h_sigma = input.initial.h_sigma;
    arma::mat h = input.initial.h;
    arma::vec h_init = input.initial.h_init;
    arma::mat h_lag = arma::zeros<arma::mat>(tt, k);

    const arma::vec h_sigma_post_shape = input.u_sigma_prior.state.sigma.shape + tt * 0.5;
    const arma::vec &h_sigma_prior_rate = input.u_sigma_prior.state.sigma.rate;
    arma::vec h_sigma_post_scale;

    const arma::vec &h0_prior_mu = input.u_sigma_prior.state.initial_state.mu;
    const arma::mat &h0_prior_v_inv = input.u_sigma_prior.state.initial_state.v_inv;
    arma::mat h0_post_v, h0_sigma_inv;

    arma::mat u_omega_inv_diag = arma::eye<arma::mat>(k * tt, k * tt);
    u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));
    arma::mat u_sigma_inv_diag = use_psi ? arma::mat(arma::trans(Psi) * u_omega_inv_diag * Psi)
                                         : u_omega_inv_diag;
    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);

    out.u_omega_inv = arma::mat(k * tt, iterations);
    out.u_sigma_inv = arma::mat(kk * tt, iterations);

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            for (int i = 0; i < tt; i++)
            {
                u_sigma.rows(k * i, k * (i + 1) - 1) = arma::solve(
                    u_sigma_inv_diag.submat(k * i, k * i, k * (i + 1) - 1, k * (i + 1) - 1),
                    diag_k);
            }

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a
            a = kalman_durbin_koopman_2002(const_cast<arma::mat &>(ymat), z, u_sigma, a_sigma,
                                           a_B, a0, a_sigma)
                    .cols(1, tt);

            // Draw a_sigma
            a_lag.col(0) = a0;
            a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
            a_lag = a - a_lag;
            a_sigma_post_scale = 1 / (a_sigma_prior_rate + arma::sum(arma::pow(a_lag, 2), 1) * 0.5);
            for (int i = 0; i < nparams; i++)
            {
                a_sigma(i, i) = 1 / arma::randg<double>(
                                        arma::distr_param(a_sigma_post_shape(i), a_sigma_post_scale(i)));
            }

            // Draw a0
            a0_sigma_inv.diag() = 1 / a_sigma.diag();
            a0_post_v = a0_prior_v_inv + a0_sigma_inv;
            a0 = draw_normal_precision(a0_post_v, a0_prior_v_inv * a0_prior_mu + a0_sigma_inv * a.col(0));

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

        // Update psi
        if (use_psi)
        {
            psi_y = arma::reshape(arma::vectorise(u.rows(1, k - 1)), k - 1, tt);
            build_psi_regressors(psi_z, u);
            for (int j = 0; j < tt; j++)
            {
                psi_u_omega.rows(j * (k - 1), (j + 1) * (k - 1) - 1) =
                    u_omega_inv_diag.submat(j * k + 1, j * k + 1, (j + 1) * k - 1, (j + 1) * k - 1);
                psi_u_omega.rows(j * (k - 1), (j + 1) * (k - 1) - 1).diag() =
                    1 / psi_u_omega.rows(j * (k - 1), (j + 1) * (k - 1) - 1).diag();
            }

            if (use_varsel_psi)
            {
                psi_z_bvs = psi_z;
                psi_z = psi_z * psi_bvs->lambda_diag;
            }

            psi = kalman_durbin_koopman_2002(psi_y, psi_z, psi_u_omega, psi_sigma, psi_B, psi0,
                                             psi_sigma)
                      .cols(1, tt);

            // Draw psi_sigma
            psi_lag.col(0) = psi0;
            psi_lag.cols(1, tt - 1) = psi.cols(0, tt - 2);
            psi_lag = psi - psi_lag;
            psi_sigma_post_scale =
                1 / (input.psi_prior.sigma.rate + arma::sum(arma::pow(psi_lag, 2), 1) * 0.5);
            for (int i = 0; i < n_psi; i++)
            {
                psi_sigma(i, i) = 1 / arma::randg<double>(arma::distr_param(
                                          psi_sigma_post_shape(i), psi_sigma_post_scale(i)));
            }

            // Draw psi0
            psi0_sigma_inv.diag() = 1 / psi_sigma.diag();
            psi0_post_v = input.psi_prior.initial_state.v_inv + psi0_sigma_inv;
            psi0 = draw_normal_precision(psi0_post_v,
                                         input.psi_prior.initial_state.v_inv * input.psi_prior.initial_state.mu + psi0_sigma_inv * psi.col(0));

            if (psi_bvs)
            {
                psi_z = psi_z_bvs;
                for (int i = 0; i < tt; i++)
                {
                    psi_u_omega_inv_diag.submat(i * (k - 1), i * (k - 1), (i + 1) * (k - 1) - 1,
                                                (i + 1) * (k - 1) - 1) =
                        u_omega_inv_diag.submat(i * k + 1, i * k + 1, (i + 1) * k - 1, (i + 1) * k - 1);
                }

                // path_row, for the reason spelled out at the same call in
                // var_tvp_gamma.cpp: element scope reached period 0 alone.
                bvs_sweep(*psi_bvs, psi, BvsScope::path_row, [&](const arma::mat &theta) {
                    for (int i = 0; i < tt; i++)
                    {
                        psi_theta_res.subvec(i * (k - 1), (i + 1) * (k - 1) - 1) =
                            psi_y.col(i) -
                            psi_z.rows(i * (k - 1), (i + 1) * (k - 1) - 1) * theta.col(i);
                    }
                    return -arma::as_scalar(arma::trans(psi_theta_res) * psi_u_omega_inv_diag *
                                            psi_theta_res) /
                           2;
                });
            }

            fill_psi_path(Psi, psi, k);
            for (int j = 0; j < tt; j++)
            {
                u.col(j) = Psi.submat(k * j, k * j, k * (j + 1) - 1, k * (j + 1) - 1) * u.col(j);
            }
        }

        // Update u_omega_inv: the ten-component mixture of Omori et al. (2007),
        // one column of log-volatility per variable.
        h = stochvol_ocsn_2007(arma::trans(u), h, h_sigma, h_init, h_y_offset);

        // Draw h_sigma
        h_lag.row(0) = arma::trans(h_init);
        h_lag.rows(1, tt - 1) = h.rows(0, tt - 2);
        h_lag = h - h_lag;
        h_sigma_post_scale =
            1 / (h_sigma_prior_rate + arma::trans(arma::sum(arma::pow(h_lag, 2))) * 0.5);
        for (int i = 0; i < k; i++)
        {
            h_sigma(i) = 1 / arma::randg<double>(
                                 arma::distr_param(h_sigma_post_shape(i), h_sigma_post_scale(i)));
        }

        // Draw h_init
        h0_sigma_inv = arma::diagmat(1 / h_sigma);
        h0_post_v = h0_prior_v_inv + h0_sigma_inv;
        h_init = draw_normal_precision(h0_post_v,
                                       h0_prior_v_inv * h0_prior_mu + h0_sigma_inv * arma::trans(h.row(0)));

        u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));

        // Update u_sigma_inv
        if (use_psi)
        {
            u_sigma_inv_diag = arma::trans(Psi) * u_omega_inv_diag * Psi;
        }
        else
        {
            u_sigma_inv_diag.diag() = u_omega_inv_diag.diag();
        }

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            if (use_a)
            {
                out.a.col(draw_pos) = arma::vectorise(a);
                out.a_sigma.col(draw_pos) = arma::vectorise(a_sigma.diag());
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_bvs->lambda;
                }
            }

            if (use_psi)
            {
                for (int i = 0; i < tt; i++)
                {
                    out.psi.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) = arma::vectorise(
                        Psi.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1));
                }

                out.psi_sigma.col(draw_pos) = arma::vectorise(psi_sigma.diag());

                if (use_varsel_psi)
                {
                    fill_strict_lower_triangle(Psi_lambda, psi_bvs->lambda);
                    out.psi_lambda.col(draw_pos) = arma::vectorise(Psi_lambda);
                }
            }

            // Measurement error
            out.u_omega_inv.col(draw_pos) = u_omega_inv_diag.diag();

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) =
                    arma::vectorise(u_sigma_inv_diag.submat(i * k, i * k, (i + 1) * k - 1,
                                                            (i + 1) * k - 1));
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarTvpStochvolSampler::forecast(const VarTvpStochvolInput &input,
                                              const VarTvpStochvolDraws &coefficients,
                                              Reporter &reporter) const
{
    const int k = input.spec.k;
    const int p = input.spec.p;
    const int h = input.spec.h;

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
    const int nparams = static_cast<int>(z.n_cols);
    const bool use_a = nparams > 0;

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("forecast regressors were supplied but posterior draws of a "
                                    "are missing");
    }

    // The caller hands over the period the forecast starts from, one column per
    // draw, as the header says.
    const arma::mat &a = coefficients.a;

    if (use_a && static_cast<int>(a.n_rows) != nparams)
    {
        throw std::invalid_argument("forecast regressors have " + std::to_string(nparams) +
                                    " columns but the posterior holds " +
                                    std::to_string(a.n_rows) + " coefficients per draw");
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
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat VarTvpStochvolSampler::log_likelihood(const VarTvpStochvolInput &input,
                                                const VarTvpStochvolDraws &coefficients) const
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

    // Every period has its own coefficients, so the regressors become block
    // diagonal and the whole path multiplies out in one go.
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        arma::sp_mat z_large(k * tt, nparams * tt);
        for (int i = 0; i < tt; i++)
        {
            z_large.submat(i * k, i * nparams, (i + 1) * k - 1, (i + 1) * nparams - 1) =
                z.rows(i * k, (i + 1) * k - 1);
        }
        u = u - z_large * coefficients.a;
    }

    arma::mat loglik(draws, tt);

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
