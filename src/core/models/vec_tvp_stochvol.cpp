// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_tvp_stochvol.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/forecast_states.h"
#include "core/models/model_support.h"
#include "core/models/noncentred_support.h"
#include "core/models/vec_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::covariance_root;
using core::pack_strict_lower_triangle;
using core::require_period_draws;
using core::require_state_mask;
using core::require_state_variances;
using core::score_vec_forecast;
using core::simulate_vec_forecast;
using core::simulates_states;
using core::step_coint_state;
using core::step_random_walk;
using core::VecForecastStep;
using core::build_psi_regressors;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::coint_state_transition;
using core::draw_coint_rho;
using core::draw_normal_precision;
using core::initial_state_variance;
using core::fill_psi_path;
using core::stacked_identity;
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

    // Which random walks are drawn non-centred: each block for itself, by
    // whether its prior names omega_v. See core/models/noncentred_support.h.
    // The cointegration space is not among them: its state variance is fixed
    // at the identity to pin beta's scale, so there is no variance to put a
    // prior on.
    const bool a_noncentred = use_a && input.a_prior.noncentred();
    const bool psi_noncentred = use_psi && input.psi_prior.noncentred();
    const bool h_noncentred = input.u_sigma_prior.state.noncentred();

    VecTvpStochvolDraws out;

    // Coefficients
    arma::mat a, a_B, a_sigma, a_lag;
    arma::vec a_sigma_post_shape, a_sigma_post_scale;
    arma::vec a0;
    arma::mat a0_post_v, a0_sigma_inv, a0_prior_v;

    // Non-centred: the signed standard deviations, the standardised path, and
    // the latest draw's ordinates at zero.
    arma::vec a_omega;
    arma::mat a_tilde;
    core::NoncentredCoefficients a_nc;

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
        a0_prior_v = initial_state_variance(input.a_prior.initial_state);
        a0_sigma_inv = a_sigma;
        a0_sigma_inv.diag() = 1 / a_sigma.diag();

        if (a_noncentred)
        {
            // The chain starts on the positive branch; the sign switch reaches
            // the other within a draw.
            a_omega = arma::sqrt(arma::vec(a_sigma.diag()));
            core::allocate_noncentred(out.a_noncentred, static_cast<arma::uword>(n_a),
                                      static_cast<arma::uword>(iterations));
        }

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
    arma::mat beta, beta_B, beta_P, beta_PtP, beta_sigma, z_b, ystar;
    arma::vec beta0;
    arma::mat beta0_post_v;

    // Not const, and neither is anything built from it: with a prior on rho the
    // state equation is rebuilt at the end of every draw.
    double rho = input.beta_prior.rho;
    const CointRhoPrior &rho_prior = input.beta_prior.rho_prior;

    if (use_beta)
    {
        beta = input.initial.beta;
        beta0 = input.initial.beta_init;
        z_b = arma::zeros<arma::mat>(k * tt, n_beta);
        ystar = ymat;

        // Fixed for the whole chain, unlike a_sigma and psi_sigma: the unit
        // state variance is what pins beta's scale against alpha's. See
        // TvpCointSpacePrior. Being the identity, it is its own inverse, which
        // is why the initial state draw below needs no inverse of it.
        beta_sigma = arma::eye<arma::mat>(n_beta, n_beta);

        // The transition is rho P, with P = I_r kron P_tau the identity unless
        // the file centres the space's marginal prior on a given one. P is fixed
        // for the chain and rho may not be, so beta_B is rebuilt from the two.
        beta_P = coint_state_transition(input.beta_prior.p_tau, rank, k_beta);
        beta_PtP = arma::trans(beta_P) * beta_P;
        beta_B = rho * beta_P;

        // beta_1 = rho P beta_0 + eta, so the prior precision picks up
        // rho^2 P'P and the right-hand side rho P' beta_1. Both are the
        // identity for the random walk.
        beta0_post_v = input.beta_prior.initial_state.v_inv + rho * rho * beta_PtP;

        out.beta = arma::mat(n_beta * tt, iterations);

        if (rho_prior.draw)
        {
            out.rho = arma::mat(1, iterations);
        }

        // The starting values of a and beta are unlikely to agree with each
        // other; the loadings' regressors belong to the beta that is actually
        // in hand before the first draw of a reads them.
        fill_z_alpha(z, beta, w_t, k, k_beta, rank, diag_k);
    }

    // Covariance block
    arma::mat psi, Psi, psi_B, psi_lag, psi_sigma, psi_u_omega, psi_y, psi_z, psi_z_bvs;
    arma::vec psi_sigma_post_shape, psi_sigma_post_scale;
    arma::vec psi0;
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
            // The measurement variance both coefficient blocks are drawn under,
            // as the stack of k x k blocks the smoother reads.
            for (int i = 0; i < tt; i++)
            {
                u_sigma.rows(k * i, k * (i + 1) - 1) = arma::solve(
                    u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1),
                    diag_k);
            }

            if (a_bvs)
            {
                z_masked = z * a_bvs->lambda_diag;
            }

            if (a_noncentred)
            {
                // The standardised path, then a0 and omega jointly, then the
                // signs; a is rebuilt from the three. The loadings' regressors
                // are this draw's beta' w, as they are for the centred draw.
                a_nc = core::draw_noncentred_path(ymat, z_a, u_sigma, u_sigma_inv_blocks,
                                                  input.a_prior, a0, a_omega, a_tilde, a);
                a_sigma.diag() = arma::square(a_omega);
            }
            else
            {
                // Update a, with a0 integrated out of the prior of the first period.
                // See initial_state_variance().
                a = kalman_durbin_koopman_2002(ymat, z_a, u_sigma, a_sigma, a_B, a0_prior_mu,
                                               a0_prior_v + a_sigma)
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
                for (int i = 0; i < n_a; i++)
                {
                    a_sigma(i, i) = 1 / arma::randg<double>(
                                            arma::distr_param(a_sigma_post_shape(i), a_sigma_post_scale(i)));
                }
            }

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

                // rho P beta0, not beta0. The smoother's sixth argument is the
                // prior mean of the state the *first* observation loads on, and
                // it does not put the transition through it -- so what belongs
                // there is beta_1's mean under beta_1 = rho P beta_0 + eta, which
                // is the state before the sample carried forward one period.
                // Passing beta_0 itself would be the random walk's answer, and
                // is only right at rho P = I; otherwise the smoother and the
                // beta_0 draw a few lines down would be fitting different
                // models, one centring beta_1 over beta_0 and the other over
                // rho P beta_0.
                beta = kalman_durbin_koopman_2002(ystar, z_b, u_sigma, beta_sigma, beta_B,
                                                  rho * beta_P * beta0, beta_sigma)
                           .cols(0, tt - 1);

                // Draw beta0
                beta0 = draw_normal_precision(
                    beta0_post_v, input.beta_prior.initial_state.v_inv *
                                          input.beta_prior.initial_state.mu +
                                      rho * arma::trans(beta_P) * beta.col(0));

                // Draw rho
                //
                // Last of the cointegration block, because it is the one
                // quantity there that conditions on the whole path and on the
                // state before it. What it moves is the state equation itself,
                // so the transition and the initial state's posterior precision
                // are rebuilt from it for the next draw.
                if (rho_prior.draw)
                {
                    rho = draw_coint_rho(beta, beta0, beta_P, rho_prior);
                    beta_B = rho * beta_P;
                    beta0_post_v = input.beta_prior.initial_state.v_inv + rho * rho * beta_PtP;
                }

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

            if (use_beta)
            {
                out.beta.col(draw_pos) = arma::vectorise(beta);
                if (rho_prior.draw)
                {
                    out.rho(0, draw_pos) = rho;
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
            out.h_sigma.col(draw_pos) = h_sigma;
            if (h_noncentred)
            {
                core::store_noncentred(out.h_noncentred, static_cast<arma::uword>(draw_pos),
                                       h_omega, h_nc);
            }

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) =
                    arma::vectorise(u_sigma_inv_blocks.rows(i * k, (i + 1) * k - 1));
            }
        }
    }

    reporter.finish();
    return out;
}

namespace
{

/// How a VecTvpStochvol's states move over a horizon, in one place.
///
/// The forecast and the score need the same walk and have to take it
/// identically, or the two would describe different models from the same file.
/// Built only on the simulated path; held states convert once and are the
/// constant VEC's case.
struct VecTvpStochvolWalk
{
    const VecTvpStochvolInput &input;
    const VecTvpStochvolDraws &coefficients;
    int k;
    arma::uword n_a;
    bool use_beta;
    bool use_psi;
    arma::mat transition;
    arma::mat diag_k;
    arma::vec a_state, a_sigma, a_mask, beta_state, psi_state, psi_sigma, psi_mask, h_state,
        h_sigma;
    double rho;

    VecTvpStochvolWalk(const VecTvpStochvolInput &in, const VecTvpStochvolDraws &draws_in)
        : input(in), coefficients(draws_in), k(in.spec.k),
          n_a(static_cast<arma::uword>(in.spec.nparams_per_period_vec())),
          use_beta(in.use_beta()), use_psi(in.use_psi()),
          diag_k(arma::eye<arma::mat>(in.spec.k, in.spec.k)), rho(in.beta_prior.rho)
    {
        const arma::uword k_u = static_cast<arma::uword>(k);
        const arma::uword draws = draws_in.iterations();
        const arma::uword n_beta = static_cast<arma::uword>(in.spec.n_beta());

        if (n_a > 0)
        {
            require_period_draws(draws_in.a, n_a, draws, "a");
            require_state_variances(draws_in.a_sigma, n_a, draws, "the coefficients");
            require_state_mask(draws_in.a_lambda, n_a, draws, "the coefficients");
        }
        if (use_beta)
        {
            require_period_draws(draws_in.beta, n_beta, draws, "beta");
            if (draws_in.has_rho())
            {
                require_period_draws(draws_in.rho, 1, draws, "rho");
            }
            transition = core::coint_state_transition(in.beta_prior.p_tau, in.spec.rank,
                                                      in.spec.k_beta);
        }
        if (use_psi)
        {
            require_period_draws(draws_in.psi, k_u * k_u, draws, "Psi");
            require_state_variances(draws_in.psi_sigma,
                                    static_cast<arma::uword>(in.spec.n_psi()), draws,
                                    "the covariance block");
            require_state_mask(draws_in.psi_lambda, k_u * k_u, draws, "the covariance block");
        }
        require_period_draws(draws_in.u_omega_inv, k_u, draws, "u_omega_inv");
        require_state_variances(draws_in.h_sigma, k_u, draws, "the log-volatilities");
    }

    void operator()(const arma::uword draw, const int i, VecForecastStep &out)
    {
        if (i == 0)
        {
            if (n_a > 0)
            {
                a_state = coefficients.a.col(draw);
                a_sigma = coefficients.a_sigma.col(draw);
                if (coefficients.a_lambda.n_elem > 0)
                {
                    a_mask = coefficients.a_lambda.col(draw);
                }
            }
            if (use_beta)
            {
                beta_state = coefficients.beta.col(draw);
                rho = coefficients.has_rho() ? coefficients.rho(0, draw) : input.beta_prior.rho;
            }
            if (use_psi)
            {
                psi_state =
                    pack_strict_lower_triangle(arma::reshape(coefficients.psi.col(draw), k, k));
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

        // The coefficients, the cointegration vectors, Psi, then the
        // log-volatilities.
        if (n_a > 0)
        {
            step_random_walk(a_state, a_sigma, a_mask);
            out.period.a = a_state;
        }
        if (use_beta)
        {
            step_coint_state(beta_state, rho, transition);
            out.period.beta = beta_state;
        }
        if (use_psi)
        {
            step_random_walk(psi_state, psi_sigma, psi_mask);
        }
        step_random_walk(h_state, h_sigma, arma::vec());

        // Psi' Omega^-1 Psi, as the sampler forms it period by period, and the
        // root of its inverse from Psi and Omega directly.
        arma::mat Psi = diag_k;
        if (use_psi)
        {
            core::fill_strict_lower_triangle(Psi, psi_state);
        }
        out.period.u_sigma_inv =
            arma::vectorise(arma::trans(Psi) * arma::diagmat(arma::exp(-h_state)) * Psi);
        out.error_root = covariance_root(Psi, arma::exp(h_state));
    }
};

} // namespace

ForecastDraws VecTvpStochvolSampler::forecast(const VecTvpStochvolInput &input,
                                              const VecTvpStochvolDraws &coefficients,
                                              Reporter &reporter) const
{
    // Everything moves with time, so the forecast starts all of it from the
    // last in-sample period -- which is what the caller is expected to have
    // sliced out, one column per draw, as the header says.
    if (!simulates_states(input.spec))
    {
        // Held: this is the constant-coefficient VEC's forecast exactly --
        // rewrite each draw in the level parameterisation and let the VAR
        // simulate the path.
        //
        // Nothing is validated here that the two calls below do not already
        // check: the conversion rejects draws that do not match the spec --
        // including the last-period slicing being wrong, which shows up as the
        // wrong row count -- and the VAR's forecast rejects a `z` that does not
        // match the converted coefficients.
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

    // Simulated: every horizon the coefficients take a step of their random
    // walk, the cointegration vectors a step of their state equation, Psi and
    // the log-volatilities steps of theirs, and the level VAR and the precision
    // are rebuilt from them. See core::simulate_vec_forecast().
    VecTvpStochvolWalk walk(input, coefficients);
    return ForecastDraws{simulate_vec_forecast(input.spec, input.forecast,
                                               coefficients.iterations(), reporter, walk)};
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
    //
    // Every period is scored under its own precision as well. The stochastic
    // volatility moves it, and scoring the whole sample under the last period's,
    // as this used to, is the likelihood of a model whose volatility does not.
    const arma::mat diag_k = arma::eye(k, k);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv, z_period;
    arma::vec resid;
    const arma::uword kk = static_cast<arma::uword>(k) * k;
    const arma::uword u_stride = core::precision_stride(coefficients.u_sigma_inv, k, tt);
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

arma::mat VecTvpStochvolSampler::predictive_log_density(
    const VecTvpStochvolInput &input, const VecTvpStochvolDraws &coefficients) const
{
    if (!simulates_states(input.spec))
    {
        // Held: the constant VEC's case exactly -- one conversion and the level
        // VAR's score, which is also where the structural refusal lives.
        VarNormalWishartInput var_input;
        var_input.spec = vec_to_var_spec(input.spec);
        var_input.forecast = input.forecast;
        var_input.test = input.test;

        VecNormalWishartDraws last_period;
        last_period.a = coefficients.a;
        last_period.beta = coefficients.beta;
        last_period.u_sigma_inv = coefficients.u_sigma_inv;

        return VarNormalWishartSampler{}.predictive_log_density(
            var_input, vec_to_var_coefficients(input.spec, last_period));
    }

    // Simulated: the same walk the forecast takes, against the realised levels
    // instead of a drawn path.
    VecTvpStochvolWalk walk(input, coefficients);
    return score_vec_forecast(input.spec, input.forecast, input.test.y,
                              coefficients.iterations(), walk);
}

} // namespace bayests
