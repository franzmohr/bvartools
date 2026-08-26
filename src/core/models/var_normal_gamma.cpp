// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_gamma.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/ssvs.h"
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
using core::SsvsBlock;
using core::ssvs_sweep;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarNormalGammaDraws VarNormalGammaSampler::draw_coefficients(const VarNormalGammaInput &input,
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

    // Both SUR blocks below need a precision that is the same in every period
    // spread over the whole sample, kron(I_tt, .): tt identical blocks, and
    // hence a density of 1/tt. Dense, the coefficient block alone would be
    // k^2 tt^2 doubles rebuilt on every draw -- 72 MB for k = 6, tt = 500 --
    // against k^2 tt nonzeros.
    const arma::sp_mat diag_tt = arma::speye<arma::sp_mat>(tt, tt);

    const bool use_psi = input.use_psi();

    const bool use_ssvs = input.spec.varsel == VarSelection::ssvs;
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_ssvs || use_bvs;

    VarNormalGammaDraws out;

    // Coefficients
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_post_v;

    // Variable selection. Only one of the two ever holds a value -- the schemes
    // are mutually exclusive -- and neither holds one when the model selects
    // nothing, which is what says "these regressors are not being selected
    // over" now that it is no longer a bool guarding a pile of empty matrices.
    std::optional<SsvsBlock> a_ssvs;
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

            if (use_ssvs)
            {
                a_ssvs.emplace(input.initial.a_lambda, input.a_varsel_prior);
            }

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
    arma::mat Psi, Psi_lambda, psi_prior_vinv, psi_post_v, psi_z, dpsi_z;
    arma::mat psi_u_omega_inv;
    arma::sp_mat psi_u_omega_inv_diag;

    std::optional<SsvsBlock> psi_ssvs;
    std::optional<BvsBlock> psi_bvs;
    arma::mat psi_z_bvs;

    if (use_psi)
    {
        n_psi = k * (k - 1) / 2;
        psi_prior_mu = input.psi_prior.mu;
        psi_prior_vinv = input.psi_prior.v_inv;
        psi = input.initial.psi;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        Psi = arma::eye<arma::mat>(k, k);
        Psi_lambda = arma::eye<arma::mat>(k, k);
        out.psi = arma::mat(k * k, iterations);

        if (use_varsel)
        {
            out.psi_lambda = arma::mat(k * k, iterations);

            if (use_ssvs)
            {
                psi_ssvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            }

            if (use_bvs)
            {
                // Captured before the first draw fills psi_z, so this is the
                // zero matrix -- and the selection step below restores it as
                // such. Kept as it stands: changing what the BVS residuals are
                // measured against would move the posterior, which is a
                // modelling decision rather than part of this split.
                psi_z_bvs = psi_z;
                psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            }
        }
    }

    // Error term
    const arma::vec u_sigma_post_shape = input.u_sigma_prior.shape + tt * 0.5;
    const arma::vec &u_sigma_prior_rate = input.u_sigma_prior.rate;
    arma::mat sse;
    arma::mat u_sigma_inv = input.initial.u_sigma_inv;
    arma::mat u_omega_inv = u_sigma_inv;
    out.u_omega_inv = arma::mat(k, iterations);
    out.u_sigma_inv = arma::mat(k * k, iterations);
    arma::sp_mat u_sigma_inv_diag;
    arma::mat dz;
    arma::mat u = arma::reshape(y, k, tt);

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            u_sigma_inv_diag = arma::kron(diag_tt, arma::sp_mat(u_sigma_inv));

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a
            //
            // The precision is symmetric, so z' D = (D z)' and one sparse
            // product serves both the posterior precision and its right-hand
            // side. Keeping the sparse operand on the left is also what picks
            // Armadillo's sparse-times-dense path rather than promoting the
            // whole block diagonal back to dense.
            dz = u_sigma_inv_diag * z;
            a_post_v = a_prior_vinv + arma::trans(dz) * z;
            a = draw_normal_precision(a_post_v,
                                      a_prior_vinv * a_prior_mu + arma::trans(dz) * y);

            if (a_ssvs)
            {
                ssvs_sweep(*a_ssvs, a, a_prior_vinv);
            }

            if (a_bvs)
            {
                z = z_bvs;
                // res' kron(I_tt, S) res is sum_t u_t' S u_t, that is
                // trace(S U U') for the k x tt error matrix U. Contracting over
                // k rather than forming the (k tt) quadratic form matters here
                // and nowhere else: this closure runs twice per selected
                // coefficient per draw, so it is the hottest arithmetic in the
                // sampler.
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::mat res = arma::reshape(y - z * theta, k, tt);
                    return -arma::accu((u_sigma_inv * res) % res) / 2;
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

            // Equation i is explained by the errors above it, so the psi block
            // carries k - 1 rows per period rather than k, and its precision is
            // the corresponding corner of u_omega_inv.
            psi_u_omega_inv = u_omega_inv.submat(1, 1, k - 1, k - 1);
            psi_u_omega_inv_diag = arma::kron(diag_tt, arma::sp_mat(psi_u_omega_inv));
            dpsi_z = psi_u_omega_inv_diag * psi_z;
            psi_post_v = psi_prior_vinv + arma::trans(dpsi_z) * psi_z;
            psi = draw_normal_precision(psi_post_v,
                                        psi_prior_vinv * psi_prior_mu + arma::trans(dpsi_z) * psi_y);

            if (psi_ssvs)
            {
                ssvs_sweep(*psi_ssvs, psi, psi_prior_vinv);
            }

            if (psi_bvs)
            {
                psi_z = psi_z_bvs;
                bvs_sweep(*psi_bvs, psi, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::mat res = arma::reshape(psi_y - psi_z * theta, k - 1, tt);
                    return -arma::accu((psi_u_omega_inv * res) % res) / 2;
                });
            }

            fill_strict_lower_triangle(Psi, psi);
            if (use_varsel)
            {
                fill_strict_lower_triangle(Psi_lambda, psi_ssvs ? psi_ssvs->lambda : psi_bvs->lambda);
            }
            u = Psi * u;
        }

        // Update u_sigma_inv
        sse = u * u.t();
        for (int i = 0; i < k; i++)
        {
            u_omega_inv(i, i) = arma::randg<double>(arma::distr_param(u_sigma_post_shape(i), 1 / arma::as_scalar(u_sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }
        if (use_psi)
        {
            u_sigma_inv = arma::trans(Psi) * u_omega_inv * Psi;
        }
        else
        {
            u_sigma_inv = u_omega_inv;
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
                    out.a_lambda.col(draw_pos) = a_ssvs ? a_ssvs->lambda : a_bvs->lambda;
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
            out.u_omega_inv.col(draw_pos) = u_omega_inv.diag();
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarNormalGammaSampler::forecast(const VarNormalGammaInput &input,
                                              const VarNormalGammaDraws &coefficients,
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

    const int nparams = (have_z || structural) ? static_cast<int>(coefficients.a.n_rows) : 0;
    const bool use_a = nparams > 0 && nparams > n_structural;

    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

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

arma::mat VarNormalGammaSampler::log_likelihood(const VarNormalGammaInput &input,
                                                const VarNormalGammaDraws &coefficients) const
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
