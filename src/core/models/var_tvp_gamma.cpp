// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_tvp_gamma.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
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

VarTvpGammaDraws VarTvpGammaSampler::draw_coefficients(const VarTvpGammaInput &input,
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
    arma::mat ymat = arma::reshape(y, k, tt);

    const bool use_psi = input.use_psi();

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VarTvpGammaDraws out;

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

    // Residuals of a candidate coefficient path, k x tt and preallocated: the
    // selection step evaluates two of them per position it revisits.
    arma::mat a_theta_res;

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
                a_theta_res = arma::zeros<arma::mat>(k, tt);
            }
        }
    }

    // Covariance block, itself time-varying
    int n_psi = 0;
    arma::mat psi, Psi, psi_B, psi_lag, psi_sigma, psi_y, psi_z, psi_z_bvs;
    arma::vec psi_sigma_post_shape, psi_sigma_post_scale;
    const arma::vec &psi_sigma_prior_rate = input.psi_prior.sigma.rate;

    arma::vec psi0;
    arma::mat psi0_post_v, psi0_sigma_inv;
    const arma::vec &psi0_prior_mu = input.psi_prior.initial_state.mu;
    const arma::mat &psi0_prior_v_inv = input.psi_prior.initial_state.v_inv;

    bool use_varsel_psi = false;
    std::optional<BvsBlock> psi_bvs;
    arma::mat psi_theta_res;
    arma::mat Psi_lambda, psi_u_omega_inv;

    if (use_psi)
    {
        n_psi = k * (k - 1) / 2;
        psi = input.initial.psi;
        psi_lag = psi;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        out.psi = arma::mat(k * k * tt, iterations);
        out.psi_sigma = arma::mat(n_psi, iterations);
        psi_B = arma::eye<arma::mat>(n_psi, n_psi);
        Psi = arma::eye<arma::mat>(k * tt, k * tt);
        fill_psi_path(Psi, psi, k);

        psi_sigma = input.initial.psi_sigma_inv;
        psi_sigma.diag() = 1 / psi_sigma.diag();
        psi_sigma_post_shape = input.psi_prior.sigma.shape + tt * 0.5;
        psi_sigma_post_scale = psi_sigma_prior_rate;

        psi0 = input.initial.psi_init;
        psi0_sigma_inv = psi_sigma;
        psi0_sigma_inv.diag() = 1 / psi_sigma.diag();

        use_varsel_psi = input.uses_psi_varsel();

        if (use_varsel_psi)
        {
            out.psi_lambda = arma::mat(k * k, iterations);
            Psi_lambda = arma::eye<arma::mat>(k, k);
            psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            psi_theta_res = arma::zeros<arma::mat>(k - 1, tt);
        }
    }

    // Error term
    const arma::vec u_sigma_post_shape = input.u_sigma_prior.shape + tt * 0.5;
    const arma::vec &u_sigma_prior_rate = input.u_sigma_prior.rate;

    arma::mat sse;
    arma::mat u = arma::reshape(y, k, tt);
    arma::mat u_omega_inv = input.initial.u_omega_inv;

    // Inverted once, before the chain starts, and reused for every psi draw
    // thereafter -- the precision is redrawn each iteration but this is not.
    // Left alone here: refreshing it would move the posterior.
    const arma::mat u_omega = arma::solve(u_omega_inv, diag_k);

    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);

    // The error precision, held as one k x k block per period rather than as the
    // (k tt) x (k tt) block diagonal this used to spell out. Every reader of it
    // is per-period already -- the smoother inverts one block at a time, the
    // selection step sums a quadratic form period by period, and the output
    // stores one block per period -- so the off-diagonal zeros were never read.
    //
    // With a covariance block they were also expensive to produce: Psi is
    // itself k tt square, so kron(I_tt, u_omega_inv) followed by
    // Psi' Omega Psi ran two dense products of order (k tt)^3 to fill a matrix
    // whose tt diagonal blocks are each just Psi_j' u_omega_inv Psi_j. That is a
    // factor tt^2 more arithmetic than the tt small products below, and 72 MB
    // of it for k = 6, tt = 500. Blocks are stacked row-wise, matching u_sigma
    // above, which the smoother reads the same way.
    arma::mat u_sigma_inv_blocks(k * tt, k);
    const auto refresh_u_sigma_inv_blocks = [&]()
    {
        for (int i = 0; i < tt; i++)
        {
            if (use_psi)
            {
                u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) =
                    arma::trans(Psi.submat(k * i, k * i, k * (i + 1) - 1, k * (i + 1) - 1)) *
                    u_omega_inv *
                    Psi.submat(k * i, k * i, k * (i + 1) - 1, k * (i + 1) - 1);
            }
            else
            {
                u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) = u_omega_inv;
            }
        }
    };
    refresh_u_sigma_inv_blocks();

    out.u_omega_inv = arma::mat(k, iterations);
    if (use_psi)
    {
        out.u_sigma_inv = arma::mat(kk * tt, iterations);
    }
    else
    {
        out.u_sigma_inv = arma::mat(kk, iterations);
    }

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            for (int i = 0; i < tt; i++)
            {
                u_sigma.rows(k * i, k * (i + 1) - 1) = arma::solve(u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1), diag_k);
            }

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a
            a = kalman_durbin_koopman_2002(ymat, z, u_sigma, a_sigma, a_B, a0, a_sigma).cols(1, tt);

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

                // The coefficients are a path, so a regressor is in or out for
                // the whole sample: switching a position off zeroes its row in
                // every period.
                //
                // The precision differs from period to period, so the quadratic
                // form is summed block by block -- sum_t r_t' S_t r_t -- rather
                // than as one dot product against a matrix that is all zeros
                // outside those blocks.
                bvs_sweep(*a_bvs, a, BvsScope::path_row, [&](const arma::mat &theta) {
                    double quadratic_form = 0.0;
                    for (int i = 0; i < tt; i++)
                    {
                        a_theta_res.col(i) =
                            ymat.col(i) - z.rows(i * k, (i + 1) * k - 1) * theta.col(i);
                        quadratic_form += arma::dot(a_theta_res.col(i),
                                                    u_sigma_inv_blocks.rows(k * i, k * (i + 1) - 1) *
                                                        a_theta_res.col(i));
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
            // Draw psi
            psi_y = arma::reshape(arma::vectorise(u.rows(1, k - 1)), k - 1, tt);
            build_psi_regressors(psi_z, u);

            if (use_varsel_psi)
            {
                psi_z_bvs = psi_z;
                psi_z = psi_z * psi_bvs->lambda_diag;
            }

            arma::mat psi_sigma_u = u_omega.submat(1, 1, k - 1, k - 1);
            psi = kalman_durbin_koopman_2002(psi_y, psi_z,
                                             psi_sigma_u,
                                             psi_sigma, psi_B, psi0, psi_sigma)
                      .cols(1, tt);

            // Draw psi_sigma
            psi_lag.col(0) = psi0;
            psi_lag.cols(1, tt - 1) = psi.cols(0, tt - 2);
            psi_lag = psi - psi_lag;
            psi_sigma_post_scale = 1 / (psi_sigma_prior_rate + arma::sum(arma::pow(psi_lag, 2), 1) * 0.5);
            for (int i = 0; i < n_psi; i++)
            {
                psi_sigma(i, i) = 1 / arma::randg<double>(arma::distr_param(psi_sigma_post_shape(i), psi_sigma_post_scale(i)));
            }

            // Draw psi0
            psi0_sigma_inv.diag() = 1 / psi_sigma.diag();
            psi0_post_v = psi0_prior_v_inv + psi0_sigma_inv;
            psi0 = draw_normal_precision(psi0_post_v,
                                         psi0_prior_v_inv * psi0_prior_mu + psi0_sigma_inv * psi.col(0));

            if (psi_bvs)
            {
                psi_z = psi_z_bvs;

                // Equation i is explained by the errors above it, so the psi
                // block carries k - 1 rows per period and its precision is the
                // corresponding corner of u_omega_inv -- the same corner in
                // every period, unlike the coefficient block above, so the
                // quadratic form collapses to trace(S R R').
                psi_u_omega_inv = u_omega_inv.submat(1, 1, k - 1, k - 1);

                // path_row, as for `a` above: `psi` is a path, n_psi x tt, and
                // excluding a contemporaneous coefficient has to exclude it in
                // every period.
                //
                // This used to pass element scope, which on a matrix is a single
                // linear index -- column-major, so with the position below n_psi
                // it reached row `pos` of period 0 and left every later period
                // untouched. Masking the drawn path on the way out was never
                // affected, so the damage was confined to the candidates the
                // sweep scores: the excluded one differed from the included one
                // in a single period out of tt, the likelihood ratio between
                // them was correspondingly close to zero, and inclusion won
                // essentially every time. The psi selection this model reports
                // was therefore not selecting -- with tt=24 it kept all three
                // free coefficients in all 80 draws of the regression fixture,
                // against 26 of 240 once the candidate spans the sample.
                bvs_sweep(*psi_bvs, psi, BvsScope::path_row, [&](const arma::mat &theta) {
                    for (int i = 0; i < tt; i++)
                    {
                        psi_theta_res.col(i) =
                            psi_y.col(i) -
                            psi_z.rows(i * (k - 1), (i + 1) * (k - 1) - 1) * theta.col(i);
                    }
                    return -arma::accu((psi_u_omega_inv * psi_theta_res) % psi_theta_res) / 2;
                });
            }
            fill_psi_path(Psi, psi, k);
            for (int j = 0; j < tt; j++)
            {
                u.col(j) = Psi.submat(k * j, k * j, k * (j + 1) - 1, k * (j + 1) - 1) * u.col(j);
            }
        }

        // Update u_omega_inv
        sse = u * u.t();
        for (int i = 0; i < k; i++)
        {
            u_omega_inv(i, i) = arma::randg<double>(arma::distr_param(u_sigma_post_shape(i), 1 / arma::as_scalar(u_sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }

        // Update u_sigma_inv
        refresh_u_sigma_inv_blocks();

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

            // Psi
            if (use_psi)
            {
                for (int i = 0; i < tt; i++)
                {
                    out.psi.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) = arma::vectorise(Psi.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1));
                }

                out.psi_sigma.col(draw_pos) = arma::vectorise(psi_sigma.diag());

                if (use_varsel_psi)
                {
                    fill_strict_lower_triangle(Psi_lambda, psi_bvs->lambda);
                    out.psi_lambda.col(draw_pos) = arma::vectorise(Psi_lambda);
                }
            }

            // Measurement error
            out.u_omega_inv.col(draw_pos) = arma::vectorise(u_omega_inv.diag());

            if (use_psi)
            {
                for (int i = 0; i < tt; i++)
                {
                    out.u_sigma_inv.submat(i * kk, draw_pos, (i + 1) * kk - 1, draw_pos) = arma::vectorise(u_sigma_inv_blocks.rows(i * k, (i + 1) * k - 1));
                }
            }
            else
            {
                out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_omega_inv);
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarTvpGammaSampler::forecast(const VarTvpGammaInput &input,
                                           const VarTvpGammaDraws &coefficients,
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
    // period of the posterior, and that has to be sliced out by size before
    // there is anything to count.
    const bool use_a = n_non_structural > 0;

    if (nparams > 0 && !coefficients.has_a())
    {
        throw std::invalid_argument("forecasting this model needs posterior draws of a, which "
                                    "are missing");
    }

    arma::mat a = coefficients.a;
    arma::mat a0;
    if (structural)
    {
        a0 = a.rows(nparams - n_structural, nparams - 1);
        if (use_a)
        {
            a = a.rows(0, nparams - n_structural - 1);
        }
    }

    const arma::uword draws = coefficients.iterations();
    const bool p_larger_than_0 = p > 0;

    arma::mat fcst = arma::zeros<arma::mat>(h * k, draws);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::vec eigval;
    arma::mat eigvec;
    arma::mat A0_inv;

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

            if (structural)
            {
                A0_inv = arma::eye<arma::mat>(k, k);
                for (int j = 1; j < k; j++)
                {
                    A0_inv.submat(j, 0, j, j - 1) = arma::trans(a0.submat(j * (j - 1) / 2, draw, (j + 1) * j / 2 - 1, draw));
                }
                A0_inv = arma::solve(A0_inv, diag_k);
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = A0_inv * fcst.submat(i * k, draw, (i + 1) * k - 1, draw);
            }
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat VarTvpGammaSampler::log_likelihood(const VarTvpGammaInput &input,
                                             const VarTvpGammaDraws &coefficients) const
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
