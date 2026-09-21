// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_stochvol.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/forecast_states.h"
#include "core/models/model_support.h"
#include "core/models/noncentred_support.h"
#include "core/models/predictive_score.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::build_psi_regressors;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::covariance_root;
using core::draw_normal_precision;
using core::pack_strict_lower_triangle;
using core::require_period_draws;
using core::require_state_mask;
using core::require_state_variances;
using core::simulates_states;
using core::step_random_walk;
using core::initial_state_variance;
using core::fill_psi_path;
using core::stacked_identity;
using core::fill_strict_lower_triangle;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarTvpStochvolDraws VarTvpStochvolSampler::draw_coefficients(const VarTvpStochvolInput &input,
                                                             Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int kk = k * k;
    const int iterations = input.spec.iterations;
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

    // Which random walks are drawn non-centred: each block for itself, by
    // whether its prior names omega_v. See core/models/noncentred_support.h.
    const bool a_noncentred = use_a && input.a_prior.noncentred();
    const bool psi_noncentred = use_psi && input.psi_prior.noncentred();
    const bool h_noncentred = input.u_sigma_prior.state.noncentred();

    VarTvpStochvolDraws out;

    // Coefficients
    arma::mat a, a_B, a_sigma, a_lag, post_a_v;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    arma::vec a0, a0_post_mu;
    arma::mat a0_post_v, a0_sigma_inv, a0_prior_v;

    // Non-centred: the signed standard deviations, the standardised path, and
    // the latest draw's ordinates at zero.
    arma::vec a_omega;
    arma::mat a_tilde;
    core::NoncentredCoefficients a_nc;

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
        a0_prior_v = initial_state_variance(input.a_prior.initial_state);
        a0_sigma_inv = a_sigma;
        a0_sigma_inv.diag() = 1 / a_sigma.diag();

        if (a_noncentred)
        {
            // The chain starts on the positive branch; the sign switch reaches
            // the other within a draw.
            a_omega = arma::sqrt(arma::vec(a_sigma.diag()));
            core::allocate_noncentred(out.a_noncentred, static_cast<arma::uword>(nparams),
                                      static_cast<arma::uword>(iterations));
        }

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
    arma::mat psi0_post_v, psi0_sigma_inv, psi0_prior_v;

    arma::vec psi_omega;
    arma::mat psi_tilde, psi_u_omega_inv;
    core::NoncentredCoefficients psi_nc;

    std::optional<BvsBlock> psi_bvs;
    arma::vec psi_theta_res;
    arma::mat Psi_lambda;

    if (use_psi)
    {
        psi = input.initial.psi;
        psi_lag = psi;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        psi_B = arma::eye<arma::mat>(n_psi, n_psi);
        Psi = stacked_identity(k, tt);
        psi_u_omega = arma::zeros<arma::mat>((k - 1) * tt, k - 1);

        out.psi = arma::mat(kk * tt, iterations);
        out.psi_sigma = arma::mat(n_psi, iterations);

        fill_psi_path(Psi, psi, k);

        psi_sigma = input.initial.psi_sigma_inv;
        psi_sigma.diag() = 1 / psi_sigma.diag();
        psi_sigma_post_shape = input.psi_prior.sigma.shape + tt * 0.5;
        psi_sigma_post_scale = input.psi_prior.sigma.rate;

        psi0 = input.initial.psi_init;
        psi0_prior_v = initial_state_variance(input.psi_prior.initial_state);
        psi0_sigma_inv = psi_sigma;
        psi0_sigma_inv.diag() = 1 / psi_sigma.diag();

        if (psi_noncentred)
        {
            psi_omega = arma::sqrt(arma::vec(psi_sigma.diag()));
            psi_u_omega_inv = arma::zeros<arma::mat>((k - 1) * tt, k - 1);
            core::allocate_noncentred(out.psi_noncentred, static_cast<arma::uword>(n_psi),
                                      static_cast<arma::uword>(iterations));
        }

        if (use_varsel_psi)
        {
            out.psi_lambda = arma::mat(kk, iterations);
            Psi_lambda = arma::eye<arma::mat>(k, k);
            psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            psi_theta_res = arma::zeros<arma::vec>((k - 1) * tt);
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

    arma::vec h_omega;
    arma::mat h_tilde;
    core::NoncentredCoefficients h_nc;
    if (h_noncentred)
    {
        h_omega = arma::sqrt(h_sigma);
        core::allocate_noncentred(out.h_noncentred, static_cast<arma::uword>(k),
                                  static_cast<arma::uword>(iterations));
    }

    // Two per-period precisions: that of each orthogonalised error, tt x k like
    // h, and the error precision itself, one k x k block per period stacked
    // row-wise. See var_tvp_gamma.cpp for why neither
    // is the (k tt) square block diagonal: every reader is per-period already,
    // and with a covariance block Psi' Omega Psi over the whole diagonal was a
    // dense product of order (k tt)^3 on every draw.
    arma::mat u_omega_inv = 1 / arma::exp(h);
    arma::mat u_sigma_inv_blocks(k * tt, k);
    const auto refresh_u_sigma_inv_blocks = [&]()
    {
        for (int i = 0; i < tt; i++)
        {
            const arma::mat omega_inv_i = arma::diagmat(u_omega_inv.row(i));
            if (use_psi)
            {
                u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) =
                    arma::trans(Psi.rows(k * i, k * (i + 1) - 1)) * omega_inv_i *
                    Psi.rows(k * i, k * (i + 1) - 1);
            }
            else
            {
                u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) = omega_inv_i;
            }
        }
    };
    refresh_u_sigma_inv_blocks();
    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);

    out.u_omega_inv = arma::mat(k * tt, iterations);
    out.u_sigma_inv = arma::mat(kk * tt, iterations);
    out.h_sigma = arma::mat(k, iterations);

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
                    u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1),
                    diag_k);
            }

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            if (a_noncentred)
            {
                // The standardised path, then a0 and omega jointly, then the
                // signs; a is rebuilt from the three.
                a_nc = core::draw_noncentred_path(ymat, z, u_sigma, u_sigma_inv_blocks,
                                                  input.a_prior, a0, a_omega, a_tilde, a);
                a_sigma.diag() = arma::square(a_omega);
            }
            else
            {
                // Update a, with a0 integrated out of the prior of the first period.
                // See initial_state_variance().
                a = kalman_durbin_koopman_2002(const_cast<arma::mat &>(ymat), z, u_sigma, a_sigma,
                                               a_B, a0_prior_mu, a0_prior_v + a_sigma)
                        .cols(0, tt - 1);

                // Draw a0, given the path it was integrated out of and before a_sigma
                // conditions on it
                a0_sigma_inv.diag() = 1 / a_sigma.diag();
                a0_post_v = a0_prior_v_inv + a0_sigma_inv;
                a0 = draw_normal_precision(a0_post_v, a0_prior_v_inv * a0_prior_mu + a0_sigma_inv * a.col(0));

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
            }

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
                    // sum_t r_t' S_t r_t, block by block.
                    double quadratic_form = 0.0;
                    for (int i = 0; i < tt; i++)
                    {
                        const arma::vec r = a_theta_res.subvec(i * k, (i + 1) * k - 1);
                        quadratic_form +=
                            arma::dot(r, u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) * r);
                    }
                    return -quadratic_form / 2;
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
                    arma::diagmat(1 / u_omega_inv.submat(j, 1, j, k - 1));
            }

            if (use_varsel_psi)
            {
                psi_z_bvs = psi_z;
                psi_z = psi_z * psi_bvs->lambda_diag;
            }

            if (psi_noncentred)
            {
                for (int j = 0; j < tt; j++)
                {
                    psi_u_omega_inv.rows(j * (k - 1), (j + 1) * (k - 1) - 1) =
                        arma::diagmat(u_omega_inv.submat(j, 1, j, k - 1));
                }
                psi_nc = core::draw_noncentred_path(psi_y, psi_z, psi_u_omega, psi_u_omega_inv,
                                                    input.psi_prior, psi0, psi_omega, psi_tilde,
                                                    psi);
                psi_sigma.diag() = arma::square(psi_omega);
            }
            else
            {
                // With psi0 integrated out of the prior of the first period, as for a
                psi = kalman_durbin_koopman_2002(psi_y, psi_z, psi_u_omega, psi_sigma, psi_B,
                                                 input.psi_prior.initial_state.mu,
                                                 psi0_prior_v + psi_sigma)
                          .cols(0, tt - 1);

                // Draw psi0, given the path and before psi_sigma conditions on it
                psi0_sigma_inv.diag() = 1 / psi_sigma.diag();
                psi0_post_v = input.psi_prior.initial_state.v_inv + psi0_sigma_inv;
                psi0 = draw_normal_precision(psi0_post_v,
                                             input.psi_prior.initial_state.v_inv * input.psi_prior.initial_state.mu + psi0_sigma_inv * psi.col(0));

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
            }

            if (psi_bvs)
            {
                psi_z = psi_z_bvs;

                // Omega is diagonal, so the quadratic form against it is a
                // weighted sum of squares: the weights are u_omega_inv's trailing
                // k - 1 columns, laid out period by period as psi_theta_res is.
                const arma::vec psi_weights =
                    arma::vectorise(arma::trans(u_omega_inv.cols(1, k - 1)));

                // path_row, for the reason spelled out at the same call in
                // var_tvp_gamma.cpp: element scope reached period 0 alone.
                bvs_sweep(*psi_bvs, psi, BvsScope::path_row, [&](const arma::mat &theta) {
                    for (int i = 0; i < tt; i++)
                    {
                        psi_theta_res.subvec(i * (k - 1), (i + 1) * (k - 1) - 1) =
                            psi_y.col(i) -
                            psi_z.rows(i * (k - 1), (i + 1) * (k - 1) - 1) * theta.col(i);
                    }
                    return -arma::dot(psi_weights, arma::square(psi_theta_res)) / 2;
                });
            }

            fill_psi_path(Psi, psi, k);
            for (int j = 0; j < tt; j++)
            {
                u.col(j) = Psi.rows(k * j, k * (j + 1) - 1) * u.col(j);
            }
        }

        if (h_noncentred)
        {
            // The same mixture, with the standardised log-volatility drawn and
            // h_init and omega regressed on it; h is rebuilt from the three.
            h_nc = core::draw_noncentred_log_volatility(arma::trans(u), input.u_sigma_prior.state,
                                                        h_y_offset, h, h_init, h_omega, h_tilde);
            h_sigma = arma::square(h_omega);
        }
        else
        {
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
        }

        u_omega_inv = 1 / arma::exp(h);

        // Update u_sigma_inv
        refresh_u_sigma_inv_blocks();

        // Store draws
        if (input.spec.keeps(draw))
        {
            const int draw_pos = input.spec.kept_index(draw);

            if (use_a)
            {
                out.a.col(draw_pos) = arma::vectorise(a);
                out.a_sigma.col(draw_pos) = arma::vectorise(a_sigma.diag());
                if (a_noncentred)
                {
                    core::store_noncentred(out.a_noncentred, static_cast<arma::uword>(draw_pos),
                                           a_omega, a_nc);
                }
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
                        Psi.rows(i * k, (i + 1) * k - 1));
                }

                out.psi_sigma.col(draw_pos) = arma::vectorise(psi_sigma.diag());
                if (psi_noncentred)
                {
                    core::store_noncentred(out.psi_noncentred, static_cast<arma::uword>(draw_pos),
                                           psi_omega, psi_nc);
                }

                if (use_varsel_psi)
                {
                    fill_strict_lower_triangle(Psi_lambda, psi_bvs->lambda);
                    out.psi_lambda.col(draw_pos) = arma::vectorise(Psi_lambda);
                }
            }

            // Measurement error
            out.u_omega_inv.col(draw_pos) = arma::vectorise(arma::trans(u_omega_inv));

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) =
                    arma::vectorise(u_sigma_inv_blocks.rows(i * k, (i + 1) * k - 1));
            }

            out.h_sigma.col(draw_pos) = h_sigma;
            if (h_noncentred)
            {
                core::store_noncentred(out.h_noncentred, static_cast<arma::uword>(draw_pos),
                                       h_omega, h_nc);
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

    arma::mat x = input.forecast.x;

    require_forecast_regressors(input.spec, x);
    core::require_forecast_horizons(x, h);

    // Counted off the model's dimensions rather than off `x`: the coefficients
    // move with time, so what the forecast starts from is the last in-sample
    // period of the posterior, and the contemporaneous block at the end of it
    // has no column in `x` to be counted by.
    const bool use_a = n_non_structural > 0;

    if (nparams > 0 && !coefficients.has_a())
    {
        throw std::invalid_argument("forecast regressors were supplied but posterior draws of a "
                                    "are missing");
    }

    // The caller hands over the period the forecast starts from, one column per
    // draw, as the header says.
    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

    if (use_a && x.n_cols * static_cast<arma::uword>(k) != a.n_rows)
    {
        throw std::invalid_argument(
            "forecast regressors and coefficient draws disagree: x has " +
            std::to_string(x.n_cols) + " columns, which over k = " + std::to_string(k) +
            " equations is " + std::to_string(x.n_cols * static_cast<arma::uword>(k)) +
            " coefficients, and a has " + std::to_string(a.n_rows) +
            " rows after the structural split");
    }

    const arma::uword draws = coefficients.iterations();
    const bool p_larger_than_0 = p > 0;

    // Whether each draw's random walks are carried over the horizon or held at
    // the end of the sample: see core/models/forecast_states.h. Simulating reads
    // more of the posterior than holding does -- how far each walk moves per
    // period, and the two halves the precision is rebuilt from at every horizon
    // -- so it is checked for up front rather than inside the loop.
    const bool simulate = simulates_states(input.spec);
    const bool use_psi = input.use_psi();
    const arma::uword k_u = static_cast<arma::uword>(k);
    if (simulate)
    {
        if (nparams > 0)
        {
            require_state_variances(coefficients.a_sigma, static_cast<arma::uword>(nparams),
                                    draws, "the coefficients");
            require_state_mask(coefficients.a_lambda, static_cast<arma::uword>(nparams), draws,
                               "the coefficients");
        }
        if (use_psi)
        {
            require_period_draws(coefficients.psi, k_u * k_u, draws, "Psi");
            require_state_variances(coefficients.psi_sigma,
                                    static_cast<arma::uword>(input.spec.n_psi()), draws,
                                    "the covariance block");
            require_state_mask(coefficients.psi_lambda, k_u * k_u, draws, "the covariance block");
        }
        require_period_draws(coefficients.u_omega_inv, k_u, draws, "u_omega_inv");
        require_state_variances(coefficients.h_sigma, k_u, draws, "the log-volatilities");
    }

    arma::mat fcst = arma::zeros<arma::mat>(h * k, draws);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::vec eigval;
    arma::mat eigvec;

    // What a simulated forecast carries from one horizon to the next.
    arma::vec a_state, a_sigma, a_mask, psi_state, psi_sigma, psi_mask, h_state, h_sigma;
    arma::mat error_root;

    // Calculate forecasts
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        // Once per draw: nothing in either depends on the horizon unless the
        // states are simulated, and then both are rebuilt at every one.
        arma::mat a0_inv =
            structural ? structural_inverse(a0, draw, diag_k) : arma::mat();
        // The draw's coefficients as the k x n_x matrix they are. The SUR
        // spelling this replaced had them as a vector and paid for the
        // reshape implicitly, once per horizon, by widening z instead.
        arma::mat a_draw =
            use_a ? arma::reshape(a.col(draw), k, x.n_cols) : arma::mat();

        if (simulate)
        {
            if (nparams > 0)
            {
                a_state = coefficients.a.col(draw);
                a_sigma = coefficients.a_sigma.col(draw);
                if (coefficients.a_lambda.n_elem > 0)
                {
                    a_mask = coefficients.a_lambda.col(draw);
                }
            }
            if (use_psi)
            {
                psi_state = pack_strict_lower_triangle(arma::reshape(coefficients.psi.col(draw), k, k));
                psi_sigma = coefficients.psi_sigma.col(draw);
                if (coefficients.psi_lambda.n_elem > 0)
                {
                    psi_mask = pack_strict_lower_triangle(
                        arma::reshape(coefficients.psi_lambda.col(draw), k, k));
                }
            }
            h_state = -arma::log(coefficients.u_omega_inv.col(draw));
            h_sigma = coefficients.h_sigma.col(draw);
        }
        else
        {
            // The error covariance factorised once per draw rather than once per
            // horizon: the precision is the same at every horizon, and the
            // factorisation draws nothing, so where it sits does not move a draw.
            arma::eig_sym(eigval, eigvec, arma::solve(arma::reshape(coefficients.u_sigma_inv.col(draw), k, k), diag_k));
        }

        for (int i = 0; i < h; i++)
        {
            if (simulate)
            {
                // Every state takes its step before the observation it generates:
                // coefficients (contemporaneous ones included), then Psi, then the
                // log-volatilities.
                if (nparams > 0)
                {
                    step_random_walk(a_state, a_sigma, a_mask);
                    if (use_a)
                    {
                        a_draw = arma::reshape(a_state.head(n_non_structural), k, x.n_cols);
                    }
                    if (structural)
                    {
                        a0_inv = structural_inverse(arma::mat(a_state.tail(n_structural)), 0, diag_k);
                    }
                }
                if (use_psi)
                {
                    step_random_walk(psi_state, psi_sigma, psi_mask);
                }
                step_random_walk(h_state, h_sigma, arma::vec());

                // Psi' Omega^-1 Psi, as the sampler forms it period by period.
                arma::mat u_sigma_inv = arma::diagmat(arma::exp(-h_state));
                if (use_psi)
                {
                    arma::mat Psi = diag_k;
                    fill_strict_lower_triangle(Psi, psi_state);
                    u_sigma_inv = arma::trans(Psi) * u_sigma_inv * Psi;
                }
                error_root = covariance_root(u_sigma_inv);
            }

            if (use_a)
            {
                // Update the lagged-endogenous columns
                if (i > 0 && p_larger_than_0)
                {
                    update_forecast_lags(x, fcst, draw, i, k, p);
                }
                // Update forecast
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = a_draw * arma::trans(x.row(i));
            }

            // Add error
            if (simulate)
            {
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) =
                    fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + error_root * arma::randn(k);
            }
            else
            {
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::trans(eigvec) * arma::randn(k);
            }

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

    // Every period under its own precision where the model's moves, and under
    // the draw's single matrix where it does not; the height of what
    // read_loglik_coefficients() hands over says which. Scoring the whole sample
    // under the last period's, as this used to, is the likelihood of a model
    // whose error covariance does not move.
    const arma::uword kk = static_cast<arma::uword>(k) * k;
    const arma::uword u_stride = core::precision_stride(coefficients.u_sigma_inv, k, tt);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    double part_b = 0.0;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (int i = 0; i < tt; i++)
        {
            if (i == 0 || u_stride != 0)
            {
                const arma::uword first = static_cast<arma::uword>(i) * u_stride;
                u_sigma_inv = arma::reshape(
                    coefficients.u_sigma_inv.submat(first, draw, first + kk - 1, draw), k, k);
                part_b = core::half_log_det_precision(u_sigma_inv);
            }
            const double part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

arma::mat VarTvpStochvolSampler::predictive_log_density(const VarTvpStochvolInput &input,
                                                       const VarTvpStochvolDraws &coefficients) const
{
    core::require_scorable(input.spec, "VarTvpStochvol");
    const arma::uword periods = core::scored_horizons(input.test.y, input.spec);
    const arma::mat x =
        core::realised_regressors(input.forecast.x, input.test.y, input.spec.k, input.spec.p);

    VarTvpStochvolInput scored;
    scored.spec = input.spec;
    scored.train.y = input.test.y.head_rows(periods);
    scored.train.x = x;
    scored.train.z = core::sur_regressors(x, input.spec.k);

    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const bool simulate = core::simulates_states(input.spec);

    VarTvpStochvolDraws scored_draws;
    scored_draws.psi = coefficients.psi;
    if (coefficients.has_a())
    {
        scored_draws.a = core::carry_state_forward(
            coefficients.a, coefficients.a_sigma, coefficients.a_lambda, periods, simulate,
            "the coefficients");
    }

    // Everything this model has moves: the coefficients above, the
    // log-volatilities, and Psi with them where the covariance block is on.
    if (simulate)
    {
        const arma::mat omega_path = core::carry_log_volatility_forward(
            coefficients.u_omega_inv, coefficients.h_sigma, periods, true);
        const arma::mat psi_path =
            input.use_psi() ? core::carry_psi_forward(coefficients.psi, coefficients.psi_sigma,
                                                      coefficients.psi_lambda, periods, k, true)
                            : arma::mat();
        scored_draws.u_omega_inv = omega_path;
        scored_draws.u_sigma_inv = core::precision_path(omega_path, psi_path, k, periods);
    }
    else
    {
        scored_draws.u_sigma_inv = arma::repmat(coefficients.u_sigma_inv, periods, 1);
        if (coefficients.u_omega_inv.n_elem > 0)
        {
            scored_draws.u_omega_inv = arma::repmat(coefficients.u_omega_inv, periods, 1);
        }
    }

    return log_likelihood(scored, scored_draws);
}

} // namespace bayests
