// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_normal_stochvol.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
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
using core::fill_strict_lower_triangle;
using core::fill_z_alpha_constant;
using core::fill_z_beta_constant;
using core::reparameterise_alpha;
using core::stacked_response;

VecNormalStochvolDraws VecNormalStochvolSampler::draw_coefficients(
    const VecNormalStochvolInput &input, Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int kk = k * k;
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;
    const arma::mat w_t = arma::trans(input.train.w);

    const int n_a = static_cast<int>(z.n_cols);
    const int rank = input.spec.rank;
    const int k_beta = input.spec.k_beta;
    const int n_alpha = input.spec.n_alpha();
    const int n_beta = input.spec.n_beta();
    const bool use_a = n_a > 0;
    const bool use_beta = input.use_beta();
    const bool use_non_alpha = n_a > n_alpha;
    const int tt = static_cast<int>(y.n_elem) / k;

    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    const arma::sp_mat diag_tt = arma::speye<arma::sp_mat>(tt, tt);

    const bool use_psi = input.use_psi();
    const int n_psi = input.spec.n_psi();

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VecNormalStochvolDraws out;

    // Coefficients (non-cointegration)
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_post_v, dz;

    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;

    if (use_a)
    {
        a_prior_mu = input.a_prior.mu;
        a_prior_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        out.a = arma::mat(n_a, iterations);

        if (use_varsel)
        {
            out.a_lambda = arma::mat(n_a, iterations);

            if (use_bvs)
            {
                z_bvs = z;
                a_bvs.emplace(input.initial.a_lambda, input.varsel_prior);
            }
        }
    }

    // Coefficients (cointegration)
    arma::mat alpha, Alpha, beta_mat, Beta_mat, BB_sqrt, diag_r;
    double coint_v_inv = 0.0;
    arma::mat coint_p_tau_inv, prior_beta_vinv, post_beta_v, z_beta, dz_beta;
    arma::vec y_beta, Beta;

    if (use_beta)
    {
        diag_r = arma::eye(rank, rank);
        beta_mat = arma::reshape(input.initial.beta, k_beta, rank);
        coint_v_inv = input.beta_prior.v_inv;
        coint_p_tau_inv = input.beta_prior.p_tau_inv;
        z_beta = arma::zeros<arma::mat>(k * tt, n_beta);
        out.beta = arma::mat(n_beta, iterations);
    }

    // Covariance block, constant even though the volatility around it is not
    arma::vec psi, psi_prior_mu, psi_y;
    arma::mat Psi, Psi_lambda, psi_prior_vinv, psi_post_v, psi_z;
    arma::sp_mat Psi_block_diagonal, psi_u_omega_inv_diag;

    std::optional<BvsBlock> psi_bvs;
    arma::mat psi_z_bvs;

    Psi = arma::eye<arma::mat>(k, k);

    if (use_psi)
    {
        psi = input.initial.psi;
        psi_prior_mu = input.psi_prior.mu;
        psi_prior_vinv = input.psi_prior.v_inv;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        Psi_lambda = arma::eye<arma::mat>(k, k);
        psi_u_omega_inv_diag = arma::speye<arma::sp_mat>((k - 1) * tt, (k - 1) * tt);
        fill_strict_lower_triangle(Psi, psi);
        out.psi = arma::mat(kk, iterations);

        if (use_varsel)
        {
            out.psi_lambda = arma::mat(kk, iterations);

            if (use_bvs)
            {
                // Captured before the first draw fills psi_z, so this is the
                // zero matrix, exactly as in var_normal_stochvol.cpp.
                psi_z_bvs = psi_z;
                psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            }
        }
    }

    // Error term
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

    arma::sp_mat u_omega_inv_diag = arma::speye<arma::sp_mat>(k * tt, k * tt);
    u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));
    arma::sp_mat u_sigma_inv_diag;

    // The precision the cointegration space prior conditions alpha on. There is
    // a different one in every period here, so the prior takes their average
    // over the sample -- the `g_i` of bvartools' .bvecalg. It appears in exactly
    // two places, alpha's prior precision and beta's, and nowhere else: the data
    // terms of both blocks use the per-period precisions in full.
    arma::mat g_i;

    // Called once before the loop as well as after every volatility draw, so the
    // chain starts from Psi' Omega^-1 Psi with the Psi the file supplied.
    // var_normal_stochvol.cpp starts from the identity instead, leaving
    // /initial/psi unread until the psi block overwrites it; deriving it here
    // makes that datum mean what it says, and moves the first draw of `a` only.
    // Called once before the loop as well as after every volatility draw, so the
    // chain starts from Psi' Omega^-1 Psi with Psi unpacked above from the
    // k(k-1)/2 free elements /initial/psi carries -- the file stores the vector,
    // never the matrix. var_normal_stochvol.cpp starts from the identity
    // instead, leaving that vector unread until the psi block overwrites it;
    // deriving it here makes the datum mean what it says, and moves the first
    // draw of `a` only.
    const auto refresh_precision = [&]()
    {
        if (use_psi)
        {
            Psi_block_diagonal = arma::kron(diag_tt, arma::sp_mat(Psi));
            u_sigma_inv_diag =
                arma::trans(Psi_block_diagonal) * u_omega_inv_diag * Psi_block_diagonal;
        }
        else
        {
            u_sigma_inv_diag = u_omega_inv_diag;
        }

        if (use_beta)
        {
            g_i = arma::mat(u_sigma_inv_diag.submat(0, 0, k - 1, k - 1));
            for (int i = 1; i < tt; i++)
            {
                g_i += arma::mat(
                    u_sigma_inv_diag.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1));
            }
            g_i /= tt;
        }
    };
    refresh_precision();

    out.u_omega_inv = arma::mat(k * tt, iterations);
    out.u_sigma_inv = arma::mat(kk * tt, iterations);

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            // Block 1: Draw non-cointegration coefficients ----
            if (use_beta)
            {
                // v^-1 (beta' P_tau^-1 beta) kron G, with G the averaged
                // precision above. Otherwise the cointegration space prior of
                // Koop, Leon-Gonzalez and Strachan (2010) exactly as
                // vec_normal_wishart.cpp writes it.
                a_prior_vinv.submat(0, 0, n_alpha - 1, n_alpha - 1) =
                    arma::kron(coint_v_inv * (arma::trans(beta_mat) * coint_p_tau_inv * beta_mat),
                               g_i);

                fill_z_alpha_constant(a_bvs ? z_bvs : z, beta_mat, w_t, n_alpha, diag_k);
            }

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a. Keeping the sparse operand on the left is what picks
            // Armadillo's sparse-times-dense path rather than promoting the
            // block diagonal to dense.
            dz = u_sigma_inv_diag * z;
            a_post_v = a_prior_vinv + arma::trans(dz) * z;
            a = draw_normal_precision(a_post_v, a_prior_vinv * a_prior_mu + arma::trans(dz) * y);

            if (a_bvs)
            {
                z = z_bvs;
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::vec res = y - z * theta;
                    return -arma::as_scalar(arma::trans(res) * u_sigma_inv_diag * res) / 2;
                });
            }
        }

        // Block 2: Draw cointegration coefficients ----
        if (use_beta)
        {
            y_beta = use_non_alpha
                         ? arma::vec(y - z.cols(n_alpha, n_a - 1) * a.subvec(n_alpha, n_a - 1))
                         : y;

            // Reparameterise alpha
            alpha = arma::reshape(a.subvec(0, n_alpha - 1), k, rank);
            Alpha = reparameterise_alpha(alpha, diag_r);

            // Update beta. Unlike vec_normal_wishart.cpp the data term cannot
            // collapse to kron(Alpha' S Alpha, sum_t w_t w_t'): that identity
            // needs one S for the whole sample, and there are tt of them here.
            // The regressors are therefore built out in full and contracted
            // against the block diagonal, which is what .bvecalg does too.
            fill_z_beta_constant(z_beta, Alpha, w_t);
            dz_beta = u_sigma_inv_diag * z_beta;

            prior_beta_vinv = arma::kron(Alpha.t() * g_i * Alpha, coint_v_inv * coint_p_tau_inv);
            post_beta_v = prior_beta_vinv + arma::trans(dz_beta) * z_beta;
            Beta = draw_normal_precision(post_beta_v, arma::trans(dz_beta) * y_beta);
            Beta_mat = arma::reshape(Beta, k_beta, rank);

            // Final cointegration values. Only the product alpha beta' is
            // identified, so the draw is split between the two by the
            // normalisation below -- and both halves have to be carried forward
            // together.
            BB_sqrt = arma::sqrtmat_sympd(arma::trans(Beta_mat) * Beta_mat);
            beta_mat = Beta_mat * arma::solve(BB_sqrt, diag_r);

            alpha = Alpha * BB_sqrt;
            a.subvec(0, n_alpha - 1) = arma::vectorise(alpha);

            fill_z_alpha_constant(z, beta_mat, w_t, n_alpha, diag_k);
            u = arma::reshape(y - z * a, k, tt);
        }
        else
        {
            u = use_a ? arma::mat(arma::reshape(y - z * a, k, tt))
                      : arma::mat(arma::reshape(y, k, tt));
        }

        // Block 3: Draw the covariance block ----
        if (use_psi)
        {
            psi_y = arma::vectorise(u.rows(1, k - 1));
            build_psi_regressors(psi_z, u);

            // The psi block explains equations 1..k-1, so its precision is the
            // per-period volatility with the first equation's row and column
            // dropped.
            for (int j = 0; j < tt; j++)
            {
                psi_u_omega_inv_diag.submat(j * (k - 1), j * (k - 1), (j + 1) * (k - 1) - 1,
                                            (j + 1) * (k - 1) - 1) =
                    u_omega_inv_diag.submat(j * k + 1, j * k + 1, (j + 1) * k - 1,
                                            (j + 1) * k - 1);
            }

            psi_post_v = psi_prior_vinv + arma::trans(psi_z) * psi_u_omega_inv_diag * psi_z;
            psi = draw_normal_precision(psi_post_v,
                                        psi_prior_vinv * psi_prior_mu +
                                            arma::trans(psi_z) * psi_u_omega_inv_diag * psi_y);

            if (psi_bvs)
            {
                psi_z = psi_z_bvs;
                bvs_sweep(*psi_bvs, psi, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::vec res = psi_y - psi_z * theta;
                    return -arma::as_scalar(arma::trans(res) * psi_u_omega_inv_diag * res) / 2;
                });
            }

            fill_strict_lower_triangle(Psi, psi);
            if (use_varsel)
            {
                fill_strict_lower_triangle(Psi_lambda, psi_bvs->lambda);
            }
            u = Psi * u;
        }

        // Block 4: Draw the log-volatility ----
        //
        // The factored routine rather than var_normal_stochvol.cpp's inline
        // copy of the same mixture: it is the one .bvecalg calls, it is what
        // test/unit_stochvol.cpp covers, and it draws the path with a banded
        // Cholesky instead of factorising a dense tt x tt precision.
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
        h_init = draw_normal_precision(
            h0_post_v, h0_prior_v_inv * h0_prior_mu + h0_sigma_inv * arma::trans(h.row(0)));

        u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));
        refresh_precision();

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            if (use_a)
            {
                out.a.col(draw_pos) = a;
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_bvs->lambda;
                }
            }

            if (use_beta)
            {
                out.beta.col(draw_pos) = arma::vectorise(beta_mat);
            }

            if (use_psi)
            {
                out.psi.col(draw_pos) = arma::vectorise(Psi);
                if (use_varsel)
                {
                    out.psi_lambda.col(draw_pos) = arma::vectorise(Psi_lambda);
                }
            }

            out.u_omega_inv.col(draw_pos) = arma::vec(u_omega_inv_diag.diag());

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) =
                    arma::vectorise(arma::mat(u_sigma_inv_diag.submat(
                        i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1)));
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VecNormalStochvolSampler::forecast(const VecNormalStochvolInput &input,
                                                 const VecNormalStochvolDraws &coefficients,
                                                 Reporter &reporter) const
{
    // Only the precision moves, so it is the only thing the caller has to have
    // sliced to the last in-sample period. From there this is the constant VEC's
    // forecast exactly -- see VecNormalWishartSampler::forecast() for why the
    // level VAR does the simulating and what that demands of
    // `input.forecast.z`.
    VarNormalWishartInput var_input;
    var_input.spec = vec_to_var_spec(input.spec);
    var_input.forecast = input.forecast;

    VecNormalWishartDraws vec_draws;
    vec_draws.a = coefficients.a;
    vec_draws.beta = coefficients.beta;
    vec_draws.u_sigma_inv = coefficients.u_sigma_inv;

    const VarNormalWishartDraws var_coefficients =
        vec_to_var_coefficients(input.spec, vec_draws);

    return VarNormalWishartSampler{}.forecast(var_input, var_coefficients, reporter);
}

arma::mat VecNormalStochvolSampler::log_likelihood(
    const VecNormalStochvolInput &input, const VecNormalStochvolDraws &coefficients) const
{
    const int k = input.spec.k;
    const int kk = k * k;

    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (coefficients.u_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
    }

    const arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;
    const arma::mat w_t = arma::trans(input.train.w);
    const int k_beta = input.spec.k_beta;
    const int rank = input.spec.rank;
    const int n_alpha = input.spec.n_alpha();
    const bool use_a = z.n_cols > 0;
    const bool use_beta = input.use_beta();
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }
    if (use_beta && !coefficients.has_beta())
    {
        throw std::invalid_argument("the model has beta coefficients but posterior draws of beta "
                                    "are missing");
    }

    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    arma::mat loglik(draws, tt);

    // Calculate errors. The coefficients are constant, so this is one product
    // per draw -- and with a cointegration relation the regressors it multiplies
    // are rebuilt from that draw's beta first.
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        if (use_beta)
        {
            for (arma::uword draw = 0; draw < draws; draw++)
            {
                fill_z_alpha_constant(
                    z, arma::reshape(coefficients.beta.col(draw), k_beta, rank), w_t, n_alpha,
                    diag_k);
                u.col(draw) = y - z * coefficients.a.col(draw);
            }
        }
        else
        {
            u = u - z * coefficients.a;
        }
    }

    // Every period has its own precision, so unlike the constant-variance VECs
    // the determinant is recomputed inside the inner loop.
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (int i = 0; i < tt; i++)
        {
            u_sigma_inv = arma::reshape(
                coefficients.u_sigma_inv.submat(i * kk, draw, (i + 1) * kk - 1, draw), k, k);
            const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
            const double part_c =
                -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) *
                                 u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) /
                2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
