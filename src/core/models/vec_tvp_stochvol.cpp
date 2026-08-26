// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_tvp_stochvol.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/model_support.h"
#include "core/models/vec_support.h"

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
using core::fill_z_alpha;
using core::fill_z_beta;
using core::stacked_response;

VecTvpStochvolDraws VecTvpStochvolSampler::draw_coefficients(const VecTvpStochvolInput &input,
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

    // k_beta x tt: the error correction term of each period is a column, which
    // is the orientation both regressor builders read it in.
    const arma::mat w_t = arma::trans(input.train.w);

    const int n_a = static_cast<int>(z.n_cols);
    const bool use_a = n_a > 0;
    const int tt = static_cast<int>(y.n_elem) / k;

    const int rank = input.spec.rank;
    const int k_beta = input.spec.k_beta;
    const int n_alpha = input.spec.n_alpha();
    const int n_beta = input.spec.n_beta();
    const bool use_beta = input.use_beta();
    const bool use_non_alpha = n_a > n_alpha;

    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::mat ymat = arma::reshape(y, k, tt);

    const bool use_psi = input.use_psi();
    const int n_psi = input.spec.n_psi();

    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;
    const bool use_varsel_psi = use_psi && input.psi_varsel == VarSelection::bvs;

    VecTvpStochvolDraws out;

    // Coefficients
    arma::mat a, a_B, a_sigma, a_lag;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    arma::vec a0;
    arma::mat a0_post_v, a0_sigma_inv;

    // Variable selection
    std::optional<BvsBlock> a_bvs;
    arma::mat z_masked;
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
        a_B = arma::eye<arma::mat>(n_a, n_a);

        out.a = arma::mat(n_a * tt, iterations);
        out.a_sigma = arma::mat(n_a, iterations);

        a_sigma_post_shape = input.a_prior.sigma.shape + tt * 0.5;
        a_sigma_post_scale = a_sigma_prior_rate;

        a0 = input.initial.a_init;
        a0_sigma_inv = a_sigma;
        a0_sigma_inv.diag() = 1 / a_sigma.diag();

        if (use_varsel)
        {
            out.a_lambda = arma::mat(n_a, iterations);

            if (use_bvs)
            {
                a_bvs.emplace(input.initial.a_lambda, input.a_varsel_prior);
                a_theta_res = arma::zeros<arma::vec>(k * tt);
            }
        }
    }

    // The regressors the smoother is handed. `z` itself when nothing is
    // selected; with BVS the masked copy, because `z` has to stay unmasked --
    // its leading columns are rebuilt from beta every draw, and the residual is
    // taken against the full matrix with the masked coefficients.
    arma::mat &z_a = a_bvs ? z_masked : z;

    // Cointegration block
    arma::mat beta, beta_B, beta_sigma, z_b, ystar;
    arma::vec beta0;
    arma::mat beta0_post_v;
    const double rho = input.beta_prior.rho;

    if (use_beta)
    {
        beta = input.initial.beta;
        beta0 = input.initial.beta_init;
        z_b = arma::zeros<arma::mat>(k * tt, n_beta);
        ystar = ymat;

        // Fixed for the whole chain, unlike a_sigma and psi_sigma: the unit
        // state variance is what pins beta's scale against alpha's. See
        // TvpCointSpacePrior. Being the identity, it is its own inverse, which
        // is why the initial state draw below reads it directly.
        beta_sigma = arma::eye<arma::mat>(n_beta, n_beta);
        beta_B = rho * arma::eye<arma::mat>(n_beta, n_beta);

        // beta_1 = rho beta_0 + eta, so the prior precision picks up rho^2 and
        // the right-hand side one rho. Both are 1 for the random walk.
        beta0_post_v = input.beta_prior.initial_state.v_inv + rho * rho * beta_sigma;

        out.beta = arma::mat(n_beta * tt, iterations);

        // The starting values of a and beta are unlikely to agree with each
        // other; the loadings' regressors belong to the beta that is actually
        // in hand before the first draw of a reads them.
        fill_z_alpha(z, beta, w_t, k, k_beta, rank, diag_k);
    }

    // Covariance block
    arma::mat psi, Psi, psi_B, psi_lag, psi_sigma, psi_u_omega, psi_y, psi_z, psi_z_bvs;
    arma::vec psi_sigma_post_shape, psi_sigma_post_scale;
    arma::vec psi0;
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
    arma::mat u = ymat;
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
            // The measurement variance both coefficient blocks are drawn under,
            // as the stack of k x k blocks the smoother reads.
            for (int i = 0; i < tt; i++)
            {
                u_sigma.rows(k * i, k * (i + 1) - 1) = arma::solve(
                    u_sigma_inv_diag.submat(k * i, k * i, k * (i + 1) - 1, k * (i + 1) - 1),
                    diag_k);
            }

            if (a_bvs)
            {
                z_masked = z * a_bvs->lambda_diag;
            }

            // Update a
            a = kalman_durbin_koopman_2002(ymat, z_a, u_sigma, a_sigma, a_B, a0, a_sigma)
                    .cols(0, tt - 1);

            // Draw a_sigma
            a_lag.col(0) = a0;
            a_lag.cols(1, tt - 1) = a.cols(0, tt - 2);
            a_lag = a - a_lag;
            a_sigma_post_scale = 1 / (a_sigma_prior_rate + arma::sum(arma::pow(a_lag, 2), 1) * 0.5);
            for (int i = 0; i < n_a; i++)
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
                // Against the unmasked `z`, as in VarTvpStochvol: the sweep
                // masks the candidate coefficients, not the regressors it is
                // scored with. The loadings are never among the positions it may
                // touch -- validate() rejects that -- so the columns beta just
                // wrote are always in.
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

            // Update beta
            if (use_beta)
            {
                // What is left of the response once everything that is not a
                // loading has been explained. The loadings' own contribution
                // stays in, because it is what beta is regressed on.
                if (use_non_alpha)
                {
                    for (int i = 0; i < tt; i++)
                    {
                        ystar.col(i) = ymat.col(i) -
                                       z.submat(i * k, n_alpha, (i + 1) * k - 1, n_a - 1) *
                                           a.submat(n_alpha, i, n_a - 1, i);
                    }
                }

                fill_z_beta(z_b, a, w_t, k, rank);

                beta = kalman_durbin_koopman_2002(ystar, z_b, u_sigma, beta_sigma, beta_B, beta0,
                                                  beta_sigma)
                           .cols(0, tt - 1);

                // Draw beta0
                beta0 = draw_normal_precision(
                    beta0_post_v, input.beta_prior.initial_state.v_inv *
                                          input.beta_prior.initial_state.mu +
                                      rho * beta.col(0));

                // Carry the new cointegration space into the regressors, for the
                // residual below and for the next draw's a block.
                fill_z_alpha(z, beta, w_t, k, k_beta, rank, diag_k);
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
                      .cols(0, tt - 1);

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

            if (use_beta)
            {
                out.beta.col(draw_pos) = arma::vectorise(beta);
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

ForecastDraws VecTvpStochvolSampler::forecast(const VecTvpStochvolInput &input,
                                              const VecTvpStochvolDraws &coefficients,
                                              Reporter &reporter) const
{
    // Everything moves with time, so the forecast holds all of it at the last
    // in-sample period -- which is what the caller is expected to have sliced
    // out, one column per draw, as the header says. From there this is the
    // constant-coefficient VEC's forecast exactly: rewrite each draw in the
    // level parameterisation and let the VAR simulate the path.
    //
    // Nothing is validated here that the two calls below do not already check:
    // the conversion rejects draws that do not match the spec -- including the
    // last-period slicing being wrong, which shows up as the wrong row count --
    // and the VAR's forecast rejects a `z` that does not match the converted
    // coefficients.
    VarNormalWishartInput var_input;
    var_input.spec = vec_to_var_spec(input.spec);
    var_input.forecast = input.forecast;

    VecNormalWishartDraws last_period;
    last_period.a = coefficients.a;
    last_period.beta = coefficients.beta;
    last_period.u_sigma_inv = coefficients.u_sigma_inv;

    const VarNormalWishartDraws var_coefficients =
        vec_to_var_coefficients(input.spec, last_period);

    return VarNormalWishartSampler{}.forecast(var_input, var_coefficients, reporter);
}

arma::mat VecTvpStochvolSampler::log_likelihood(const VecTvpStochvolInput &input,
                                                const VecTvpStochvolDraws &coefficients) const
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
    const arma::mat w_t = arma::trans(input.train.w);

    const int n_a = static_cast<int>(z.n_cols);
    const bool use_a = n_a > 0;
    const int rank = input.spec.rank;
    const int k_beta = input.spec.k_beta;
    const int n_alpha = input.spec.n_alpha();
    const int n_beta = input.spec.n_beta();
    const bool use_beta = input.use_beta();

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }
    if (use_beta && !coefficients.has_beta())
    {
        throw std::invalid_argument("the model has a cointegration relation but posterior draws of "
                                    "beta are missing");
    }

    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    const arma::mat ymat = arma::reshape(y, k, tt);

    arma::mat loglik(draws, tt);

    // Unlike the VAR's, this residual cannot be one matrix product: the leading
    // columns of the regressors are beta' w_{t-1}, so they differ by period and
    // by draw. The block is rebuilt in place, into a copy of the period's
    // regressors, and the data columns are copied along with it.
    const arma::mat diag_k = arma::eye(k, k);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv, z_period;
    arma::vec resid;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, k);
        const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;

        for (int i = 0; i < tt; i++)
        {
            resid = ymat.col(i);

            if (use_a)
            {
                z_period = z.rows(i * k, (i + 1) * k - 1);

                if (use_beta)
                {
                    z_period.cols(0, n_alpha - 1) = arma::kron(
                        arma::trans(arma::trans(arma::reshape(
                                        coefficients.beta.submat(i * n_beta, draw,
                                                                 (i + 1) * n_beta - 1, draw),
                                        k_beta, rank)) *
                                    w_t.col(i)),
                        diag_k);
                }

                resid -= z_period * coefficients.a.submat(i * n_a, draw, (i + 1) * n_a - 1, draw);
            }

            const double part_c =
                -arma::as_scalar(arma::trans(resid) * u_sigma_inv * resid) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
