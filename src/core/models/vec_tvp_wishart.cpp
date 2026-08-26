// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_tvp_wishart.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/wishart.h"
#include "core/models/model_support.h"
#include "core/models/vec_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_normal_precision;
using core::fill_z_alpha;
using core::fill_z_beta;
using core::stacked_response;

VecTvpWishartDraws VecTvpWishartSampler::draw_coefficients(const VecTvpWishartInput &input,
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

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VecTvpWishartDraws out;

    // Coefficients: a whole path per draw, moved by the simulation smoother.
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

        // Fixed for the whole chain, unlike a_sigma: the unit state variance is
        // what pins beta's scale against alpha's. See TvpCointSpacePrior. Being
        // the identity, it is its own inverse, which is why the initial state
        // draw below reads it directly on both sides.
        beta_sigma = arma::eye<arma::mat>(n_beta, n_beta);
        beta_B = rho * arma::eye<arma::mat>(n_beta, n_beta);

        // beta_1 = rho beta_0 + eta, so the prior precision picks up rho^2 and
        // the right-hand side one rho.
        beta0_post_v = input.beta_prior.initial_state.v_inv + rho * rho * beta_sigma;

        out.beta = arma::mat(n_beta * tt, iterations);

        // The starting values of a and beta are unlikely to agree with each
        // other; the loadings' regressors belong to the beta that is actually
        // in hand before the first draw of a reads them.
        fill_z_alpha(z, beta, w_t, k, k_beta, rank, diag_k);
    }

    // Error term
    const int post_u_sigma_df = input.u_sigma_prior.df + tt;
    const arma::mat &prior_u_sigma_scale = input.u_sigma_prior.scale;
    arma::mat u_sigma_inv = input.initial.u_sigma_inv;
    arma::mat u = ymat;
    out.u_sigma_inv = arma::mat(kk, iterations);

    // The measurement variance the smoother reads, one k x k block per period.
    // Constant here, so it is one inverse repeated -- but rebuilt every draw all
    // the same, since u_sigma_inv is redrawn at the end of each iteration.
    arma::mat u_sigma = arma::zeros<arma::mat>(k * tt, k);
    arma::sp_mat u_sigma_inv_diag;
    const arma::sp_mat diag_tt_sp = arma::speye<arma::sp_mat>(tt, tt);

    // Start simulation
    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (use_a)
        {
            u_sigma = arma::repmat(arma::solve(u_sigma_inv, diag_k), tt, 1);
            u_sigma_inv_diag = arma::kron(diag_tt_sp, arma::sp_mat(u_sigma_inv));

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
                // Against the unmasked `z`, as in VarTvpWishart: the sweep masks
                // the candidate coefficients, not the regressors it is scored
                // with. The loadings are never among the positions it may touch
                // -- validate() rejects that -- so the columns beta wrote are
                // always in.
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

        // Update u_sigma_inv.
        //
        // The prior scale alone, with no cointegration term added to it and no
        // rank added to the degrees of freedom -- unlike VecNormalWishart, which
        // does both. There is nothing here to add: the constant VEC's prior
        // conditions alpha on Sigma, so Sigma's posterior owes it a term back,
        // whereas the loadings here are a random walk whose innovation variance
        // is drawn from a gamma of its own and never sees Sigma.
        u_sigma_inv = wishart(u, prior_u_sigma_scale, post_u_sigma_df);

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

            out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VecTvpWishartSampler::forecast(const VecTvpWishartInput &input,
                                             const VecTvpWishartDraws &coefficients,
                                             Reporter &reporter) const
{
    // The coefficients move and the precision does not, so the caller is
    // expected to have sliced the last in-sample period out of `a` and `beta`
    // and to pass `u_sigma_inv` whole. From there this is the constant VEC's
    // forecast exactly: rewrite each draw in the level parameterisation and let
    // the VAR simulate the path. Nothing is validated here that the two calls
    // below do not already check.
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

arma::mat VecTvpWishartSampler::log_likelihood(const VecTvpWishartInput &input,
                                               const VecTvpWishartDraws &coefficients) const
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

            const double part_c = -arma::as_scalar(arma::trans(resid) * u_sigma_inv * resid) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
