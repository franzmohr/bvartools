// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_stochvol.h"

#include "core/algorithms/bvs.h"
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
using core::fill_strict_lower_triangle;
using core::fill_strict_lower_triangle_by_column;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;

namespace
{

// The ten-component normal mixture of Omori, Chib, Shephard and Nakajima
// (2007) that approximates the log chi-squared distribution of log(u^2). It is
// what turns the non-linear measurement equation of a stochastic volatility
// model into a conditionally linear one the Kalman machinery can handle.
const arma::rowvec::fixed<10> kMixtureWeight = {0.00609, 0.04775, 0.13057, 0.20674, 0.22715,
                                                0.18842, 0.12047, 0.05591, 0.01575, 0.00115};
const arma::rowvec::fixed<10> kMixtureMean = {1.92677, 1.34744, 0.73504, 0.02266, -0.85173,
                                              -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const arma::rowvec::fixed<10> kMixtureVariance = {0.11265, 0.17788, 0.26768, 0.40611, 0.62699,
                                                  0.98583, 1.57469, 2.54498, 4.16591, 7.33342};

} // namespace

VarNormalStochvolDraws VarNormalStochvolSampler::draw_coefficients(
    const VarNormalStochvolInput &input, Reporter &reporter) const
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
    const arma::sp_mat diag_tt = arma::eye<arma::sp_mat>(tt, tt);

    const bool use_psi = input.use_psi();

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VarNormalStochvolDraws out;

    // Coefficients
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_post_v;

    // Variable selection
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;

    if (use_a)
    {
        a_prior_mu = input.a_prior.mu;
        a_prior_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        out.a = arma::mat(nparams, iterations);

        if (use_varsel)
        {
            out.a_lambda = arma::mat(nparams, iterations);

            if (use_bvs)
            {
                z_bvs = z;
                a_bvs.emplace(input.initial.a_lambda, input.a_varsel_prior);
            }
        }
    }

    // Covariance block
    int n_psi = 0;
    arma::vec psi, psi_prior_mu, psi_y;
    arma::mat Psi_lambda, psi_prior_vinv, psi_post_v, psi_z;
    arma::sp_mat Psi, Psi_block_diagonal, psi_u_omega_inv_diag;

    std::optional<BvsBlock> psi_bvs;
    arma::mat psi_z_bvs;

    if (use_psi)
    {
        n_psi = k * (k - 1) / 2;
        psi = input.initial.psi;
        psi_prior_mu = input.psi_prior.mu;
        psi_prior_vinv = input.psi_prior.v_inv;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        Psi = arma::eye<arma::mat>(k, k);
        Psi_lambda = arma::eye<arma::mat>(k, k);
        psi_u_omega_inv_diag = arma::eye<arma::sp_mat>((k - 1) * tt, (k - 1) * tt);
        out.psi = arma::mat(k * k, iterations);

        if (use_varsel)
        {
            out.psi_lambda = arma::mat(k * k, iterations);

            if (use_bvs)
            {
                // Captured before the first draw fills psi_z, so this is the
                // zero matrix -- and the selection step below restores it as
                // such. Left as it stands: what the BVS residuals are measured
                // against is a modelling decision, not part of this split.
                psi_z_bvs = psi_z;
                psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            }
        }
    }

    // Error term
    arma::mat u = arma::reshape(y, k, tt);
    arma::mat h_y;
    const arma::vec &h_y_offset = input.u_sigma_prior.offset;

    // The second difference operator of the random walk the log-volatility
    // follows, as a precision: hh = D'D with D the first difference.
    arma::sp_mat hh = arma::eye<arma::sp_mat>(tt, tt);
    hh.diag(-1) = -arma::ones<arma::vec>(tt - 1);
    hh = hh.t() * hh;

    const arma::mat p_i_matrix = arma::repmat(kMixtureWeight, tt, 1);
    const arma::mat mu_matrix = arma::repmat(kMixtureMean, tt, 1);
    const arma::mat sigma_matrix = arma::repmat(arma::sqrt(kMixtureVariance), tt, 1);
    const arma::vec vec_tt = arma::ones<arma::vec>(tt);

    arma::mat h_q, sigh_hh, h_post_v;
    arma::sp_mat sigs = arma::eye<arma::sp_mat>(tt, tt);
    arma::umat s;

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

    arma::sp_mat u_omega_inv_diag = arma::eye<arma::sp_mat>(k * tt, k * tt);
    u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));
    arma::sp_mat u_sigma_inv_diag;
    if (use_psi)
    {
        Psi_block_diagonal = arma::kron(diag_tt, Psi);
        u_sigma_inv_diag = arma::trans(Psi_block_diagonal) * u_omega_inv_diag * Psi_block_diagonal;
    }
    else
    {
        u_sigma_inv_diag = u_omega_inv_diag;
    }

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

            // Update a
            a_post_v = a_prior_vinv + arma::trans(z) * u_sigma_inv_diag * z;
            a = draw_normal_precision(a_post_v,
                                      a_prior_vinv * a_prior_mu + arma::trans(z) * u_sigma_inv_diag * y);

            if (a_bvs)
            {
                z = z_bvs;
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::vec res = y - z * theta;
                    return -arma::as_scalar(arma::trans(res) * u_sigma_inv_diag * res) / 2;
                });
            }

            u = arma::reshape(y - z * a, k, tt);
        }
        else
        {
            u = arma::reshape(y, k, tt);
        }

        // Update psi
        if (use_psi)
        {
            psi_y = arma::vectorise(u.rows(1, k - 1));
            build_psi_regressors(psi_z, u);

            // The psi block explains equations 1..k-1, so its precision is the
            // per-period volatility with the first equation's row and column
            // dropped. Written once per period: the fill used to sit inside the
            // regressor loop above, which repeated each block k-1 times.
            for (int j = 0; j < tt; j++)
            {
                psi_u_omega_inv_diag.submat(j * (k - 1),
                                            j * (k - 1),
                                            (j + 1) * (k - 1) - 1,
                                            (j + 1) * (k - 1) - 1) =
                    u_omega_inv_diag.submat(j * k + 1, j * k + 1,
                                            (j + 1) * k - 1, (j + 1) * k - 1);
            }

            psi_post_v = psi_prior_vinv + arma::trans(psi_z) * psi_u_omega_inv_diag * psi_z;
            psi = draw_normal_precision(psi_post_v,
                                        psi_prior_vinv * psi_prior_mu + arma::trans(psi_z) * psi_u_omega_inv_diag * psi_y);

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

        // Update the log-volatility, one variable at a time
        h_y = arma::trans(u);
        for (int i = 0; i < k; i++)
        {
            // Prepare series
            h_y.col(i) = arma::log(arma::pow(h_y.col(i), 2) + h_y_offset(i));

            // Sample s
            h_q = p_i_matrix % arma::normpdf(arma::repmat(h_y.col(i), 1, 10), arma::repmat(h.col(i), 1, 10) + mu_matrix, sigma_matrix);
            h_q = h_q / arma::repmat(arma::sum(h_q, 1), 1, 10);
            s = 10 - arma::sum(arma::repmat(arma::randu<arma::vec>(tt), 1, 10) < arma::cumsum(h_q, 1), 1);

            // Sample log-volatility
            sigh_hh = hh / h_sigma(i);
            sigs.diag() = 1 / kMixtureVariance.elem(s);
            h_post_v = sigh_hh + sigs;
            h.col(i) = draw_normal_precision(h_post_v,
                                             sigh_hh * vec_tt * h_init(i) + sigs * (h_y.col(i) - kMixtureMean.elem(s)));
        }

        // Draw sigma_h
        h_lag.row(0) = h_init.t();
        h_lag.rows(1, tt - 1) = h.rows(0, tt - 2);
        h_lag = h - h_lag;
        h_sigma_post_scale = 1 / (h_sigma_prior_rate + arma::trans(arma::sum(arma::pow(h_lag, 2))) * 0.5);
        for (int i = 0; i < k; i++)
        {
            h_sigma(i) = 1 / arma::randg<double>(arma::distr_param(h_sigma_post_shape(i), h_sigma_post_scale(i)));
        }

        // Draw h_init
        h0_sigma_inv = arma::diagmat(1 / h_sigma);
        h0_post_v = h0_prior_v_inv + h0_sigma_inv;
        h_init = draw_normal_precision(h0_post_v,
                                       h0_prior_v_inv * h0_prior_mu + h0_sigma_inv * h.row(0).t());

        u_omega_inv_diag.diag() = 1 / arma::exp(arma::vectorise(arma::trans(h)));

        if (use_psi)
        {
            Psi_block_diagonal = arma::kron(diag_tt, Psi);
            u_sigma_inv_diag = arma::trans(Psi_block_diagonal) * u_omega_inv_diag * Psi_block_diagonal;
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
                out.a.col(draw_pos) = a;
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_bvs->lambda;
                }
            }

            if (use_psi)
            {
                out.psi.col(draw_pos) = arma::vectorise(Psi);
                if (use_varsel)
                {
                    out.psi_lambda.col(draw_pos) = arma::vectorise(Psi_lambda);
                }
            }

            out.u_omega_inv.col(draw_pos) = u_omega_inv_diag.diag();

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * k * k, draw_pos, (i + 1) * k * k - 1, draw_pos) = arma::vectorise(u_sigma_inv_diag.submat(i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1));
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarNormalStochvolSampler::forecast(const VarNormalStochvolInput &input,
                                                 const VarNormalStochvolDraws &coefficients,
                                                 Reporter &reporter) const
{
    const int k = input.spec.k;
    const int p = input.spec.p;
    const int h = input.spec.h;
    const bool structural = input.spec.structural;
    const int n_structural = input.spec.n_structural();

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

    // The coefficient draws are only consulted when there are regressors to
    // apply them to or a contemporaneous matrix to split off; without either,
    // the path is the error process alone.
    const bool have_z = z.n_elem > 0;
    if ((have_z || structural) && !coefficients.has_a())
    {
        throw std::invalid_argument("forecasting from these regressors needs posterior draws of "
                                    "a, which are missing");
    }

    // Counted off the posterior, not off z. The structural coefficients are the
    // last n_structural rows of a and have no column in z, so z.n_cols is short
    // by exactly that many: splitting on it cuts a in the wrong place, takes the
    // contemporaneous block from the lag coefficients, and leaves a column count
    // that no longer matches z.
    const int nparams = (have_z || structural) ? static_cast<int>(coefficients.a.n_rows) : 0;
    const bool use_a = nparams > 0 && nparams > n_structural;

    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

    // The invariant the split has to preserve. Checked here so a future mismatch
    // names both sides instead of surfacing as an Armadillo dimension error from
    // inside the draw loop.
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

arma::mat VarNormalStochvolSampler::log_likelihood(const VarNormalStochvolInput &input,
                                                   const VarNormalStochvolDraws &coefficients) const
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
    const arma::mat &z = input.train.z;
    const bool use_a = z.n_cols > 0;

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }

    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    arma::mat loglik = arma::mat(draws, tt);

    // Calculate errors
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        u = u - z * coefficients.a;
    }

    // Calculate log likelihood. Every period has its own precision, so unlike
    // the constant-variance models the determinant has to be recomputed inside
    // the inner loop rather than once per draw.
    const arma::mat diag_k = arma::eye(k, k);
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (int i = 0; i < tt; i++)
        {
            u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.submat(i * kk, draw, (i + 1) * kk - 1, draw), k, k);
            const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
            const double part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
