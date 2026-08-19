// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_wishart.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/ssvs.h"
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
using core::SsvsBlock;
using core::ssvs_sweep;
using core::stacked_response;
using core::require_forecast_regressors;
using core::update_forecast_lags;

VarNormalWishartDraws VarNormalWishartSampler::draw_coefficients(const VarNormalWishartInput &input,
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

    const bool use_ssvs = input.spec.varsel == VarSelection::ssvs;
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_ssvs || use_bvs;

    VarNormalWishartDraws out;

    // Coefficients
    arma::vec a, prior_a_mu;
    arma::mat prior_a_vinv;
    arma::mat post_a_v, dz;

    // The error precision is the same in every period, so the block diagonal
    // the SUR form needs is kron(I_tt, u_sigma_inv): tt identical k x k blocks,
    // and hence a density of 1/tt. Dense, it would be k^2 tt^2 doubles rebuilt
    // on every draw -- 72 MB for k = 6, tt = 500 -- against k^2 tt nonzeros.
    arma::sp_mat diag_tt, u_sigma_inv_diag;

    // Variable selection. Only one of the two is ever engaged.
    std::optional<SsvsBlock> a_ssvs;
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;

    if (use_a)
    {
        diag_tt = arma::speye<arma::sp_mat>(tt, tt);
        prior_a_mu = input.a_prior.mu;
        prior_a_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        out.a = arma::mat(nparams, iterations);

        if (use_varsel)
        {
            out.a_lambda = arma::mat(nparams, iterations);

            if (use_ssvs)
            {
                a_ssvs.emplace(input.initial.a_lambda, input.varsel_prior);
            }

            if (use_bvs)
            {
                z_bvs = z;
                a_bvs.emplace(input.initial.a_lambda, input.varsel_prior);
            }
        }
    }

    // Error term
    const int post_u_sigma_df = input.u_sigma_prior.df + tt;
    const arma::mat &prior_u_sigma_scale = input.u_sigma_prior.scale;
    arma::mat u_sigma_inv = input.initial.u_sigma_inv;
    out.u_sigma_inv = arma::mat(k * k, iterations);
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
            post_a_v = prior_a_vinv + arma::trans(dz) * z;
            a = draw_normal_precision(post_a_v,
                                      prior_a_vinv * prior_a_mu + arma::trans(dz) * y);

            if (a_ssvs)
            {
                ssvs_sweep(*a_ssvs, a, prior_a_vinv);
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

        // Update u_sigma_inv
        u_sigma_inv = wishart(u, prior_u_sigma_scale, post_u_sigma_df);

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
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarNormalWishartSampler::forecast(const VarNormalWishartInput &input,
                                                const VarNormalWishartDraws &coefficients,
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
    if (use_a && static_cast<int>(coefficients.a.n_rows) != nparams)
    {
        throw std::invalid_argument("forecast regressors have " + std::to_string(nparams) +
                                    " columns but the posterior holds " +
                                    std::to_string(coefficients.a.n_rows) + " coefficients");
    }
    if (use_a && static_cast<int>(z.n_rows) != h * k)
    {
        throw std::invalid_argument("forecast regressors must have " + std::to_string(h * k) +
                                    " rows, got " + std::to_string(z.n_rows));
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
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = z.rows(i * k, (i + 1) * k - 1) * coefficients.a.col(draw);
            }

            // Add error
            arma::eig_sym(eigval, eigvec, arma::solve(arma::reshape(coefficients.u_sigma_inv.col(draw), k, k), diag_k));
            fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::trans(eigvec) * arma::randn(k);
        }
    }

    reporter.finish();
    return ForecastDraws{fcst};
}

arma::mat VarNormalWishartSampler::log_likelihood(const VarNormalWishartInput &input,
                                                  const VarNormalWishartDraws &coefficients) const
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
