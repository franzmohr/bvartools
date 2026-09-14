// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_stochvol.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/forecast_states.h"
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
using core::covariance_root;
using core::draw_normal_precision;
using core::require_period_draws;
using core::require_state_variances;
using core::simulates_states;
using core::step_random_walk;
using core::draw_stochvol_state;
using core::fill_strict_lower_triangle;
using core::fill_strict_lower_triangle_by_column;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarNormalStochvolDraws VarNormalStochvolSampler::draw_coefficients(
    const VarNormalStochvolInput &input, Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int iterations = input.spec.iterations;
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

    const arma::vec h_sigma_post_shape = input.u_sigma_prior.state.sigma.shape + tt * 0.5;
    const arma::vec &h_sigma_prior_rate = input.u_sigma_prior.state.sigma.rate;

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
    out.h_sigma = arma::mat(k, iterations);

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

            // BVS draws psi against the regressors masked by the indicators the
            // last sweep left, and scores its candidates against the unmasked
            // ones -- kept here, once they hold this draw's errors. See
            // var_normal_gamma.cpp.
            if (psi_bvs)
            {
                psi_z_bvs = psi_z;
                psi_z = psi_z * psi_bvs->lambda_diag;
            }

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

        // Update the log-volatility ----
        //
        // The factored routine the other stochastic volatility samplers use.
        // This model carried its own copy of the ten-component mixture, and
        // the copy had both faults stochvol_mixture.h describes: component
        // probabilities formed as densities, which underflow to a row of NaNs
        // for an observation far out in the tails, and an indicator index left
        // unclamped, so that row indexed one past the end of the table. It also
        // factorised a dense tt x tt precision per variable per draw, where the
        // routine takes a banded Cholesky.
        h = stochvol_ocsn_2007(arma::trans(u), h, h_sigma, h_init, h_y_offset);
        draw_stochvol_state(h_sigma, h_init, h, h_sigma_post_shape, h_sigma_prior_rate,
                            input.u_sigma_prior.state.initial_state);

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
        if (input.spec.keeps(draw))
        {
            const int draw_pos = input.spec.kept_index(draw);

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
            out.h_sigma.col(draw_pos) = h_sigma;

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

    arma::mat x = input.forecast.x;

    require_forecast_regressors(input.spec, x);
    core::require_forecast_horizons(x, h);

    // The coefficient draws are only consulted when there are regressors to
    // apply them to or a contemporaneous matrix to split off; without either,
    // the path is the error process alone.
    const bool have_x = x.n_elem > 0;
    if ((have_x || structural) && !coefficients.has_a())
    {
        throw std::invalid_argument("forecasting from these regressors needs posterior draws of "
                                    "a, which are missing");
    }

    // Counted off the posterior, not off x. The structural coefficients are the
    // last n_structural rows of a and have no column in x, so k * x.n_cols is
    // short by exactly that many: splitting on it cuts a in the wrong place,
    // takes the contemporaneous block from the lag coefficients, and leaves a
    // width that no longer matches x.
    const int nparams = (have_x || structural) ? static_cast<int>(coefficients.a.n_rows) : 0;
    const bool use_a = nparams > 0 && nparams > n_structural;

    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

    // The invariant the split has to preserve. Checked here so a future mismatch
    // names both sides instead of surfacing as an Armadillo dimension error from
    // inside the draw loop.
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

    // Whether each draw's log-volatilities are carried over the horizon or held
    // at the end of the sample: see core/models/forecast_states.h. The
    // coefficients and Psi are constant here, so the volatility is all that
    // moves -- and with it the precision, rebuilt at every horizon.
    const bool simulate = simulates_states(input.spec);
    const bool use_psi = input.use_psi();
    const arma::uword k_u = static_cast<arma::uword>(k);
    if (simulate)
    {
        if (use_psi)
        {
            require_period_draws(coefficients.psi, k_u * k_u, draws, "Psi");
        }
        require_period_draws(coefficients.u_omega_inv, k_u, draws, "u_omega_inv");
        require_state_variances(coefficients.h_sigma, k_u, draws, "the log-volatilities");
    }

    arma::mat fcst = arma::zeros<arma::mat>(h * k, draws);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::vec eigval;
    arma::mat eigvec;

    // What a simulated forecast carries from one horizon to the next.
    arma::vec h_state, h_sigma;
    arma::mat error_root;

    // Calculate forecasts
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        // Once per draw: nothing in either depends on the horizon.
        const arma::mat a0_inv =
            structural ? structural_inverse(a0, draw, diag_k) : arma::mat();
        // The draw's coefficients as the k x n_x matrix they are. The SUR
        // spelling this replaced had them as a vector and paid for the
        // reshape implicitly, once per horizon, by widening z instead.
        const arma::mat a_draw =
            use_a ? arma::reshape(a.col(draw), k, x.n_cols) : arma::mat();

        if (simulate)
        {
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
                // The step comes before the observation it generates.
                step_random_walk(h_state, h_sigma, arma::vec());

                // Psi' Omega^-1 Psi, as the sampler forms it.
                arma::mat u_sigma_inv = arma::diagmat(arma::exp(-h_state));
                if (use_psi)
                {
                    const arma::mat Psi = arma::reshape(coefficients.psi.col(draw), k, k);
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
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (int i = 0; i < tt; i++)
        {
            u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.submat(i * kk, draw, (i + 1) * kk - 1, draw), k, k);
            const double part_b = core::half_log_det_precision(u_sigma_inv);
            const double part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
