// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_normal_wishart.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
#include "core/algorithms/ssvs.h"
#include "core/algorithms/wishart.h"
#include "core/models/model_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

    using core::bvs_sweep;
    using core::BvsBlock;
    using core::BvsScope;
    using core::draw_normal_precision;
    using core::ssvs_sweep;
    using core::SsvsBlock;
    using core::stacked_response;

    VecNormalWishartDraws VecNormalWishartSampler::draw_coefficients(const VecNormalWishartInput &input,
                                                                     Reporter &reporter) const
    {
        input.validate();

        const int k = input.spec.k;
        const int iterations = input.spec.iterations;
        const int burnin = input.spec.burnin;
        const int draws = input.spec.draws();

        const arma::vec y = stacked_response(input.train);
        arma::mat z = input.train.z;
        arma::mat w_t = arma::trans(input.train.w);

        const int n_a = static_cast<int>(z.n_cols);
        const int rank = static_cast<int>(input.spec.rank);
        const int k_beta = static_cast<int>(input.spec.k_beta);
        const int n_alpha = static_cast<int>(input.spec.n_alpha());
        const bool use_non_alpha = n_a > n_alpha;
        const int n_beta = static_cast<int>(input.spec.n_beta());
        const bool use_a = n_a > 0;
        const bool use_beta = rank > 0;
        const int tt = static_cast<int>(y.n_elem) / k;

        const bool use_ssvs = input.spec.varsel == VarSelection::ssvs;
        const bool use_bvs = input.spec.varsel == VarSelection::bvs;
        const bool use_varsel = use_ssvs || use_bvs;

        VecNormalWishartDraws out;

        // Coefficients (non-cointegration)
        arma::vec a, prior_a_mu;
        arma::mat prior_a_vinv, diag_k;
        arma::mat post_a_v, dz;

        // The error precision is the same in every period, so the block
        // diagonal the SUR form needs is kron(I_tt, u_sigma_inv): tt identical
        // k x k blocks, and hence a density of 1/tt. Dense, it would be
        // k^2 tt^2 doubles rebuilt on every draw -- 72 MB for k = 6, tt = 500 --
        // against k^2 tt nonzeros. Same choice, for the same reason, as in
        // var_tvp_wishart.cpp.
        arma::sp_mat diag_tt, u_sigma_inv_diag;

        // Coefficients (cointegration)
        arma::mat alpha, Alpha, beta, beta_mat, Beta, Beta_mat, BB_sqrt, diag_r;
        // Initialised, not merely declared: block 3 reads it whenever the model
        // has a cointegration relation, and a rank-zero VEC never assigns it.
        double coint_v_inv = 0.0;
        arma::mat prior_beta_vinv, post_beta_v;
        arma::vec y_beta;
        arma::mat coint_p_tau_inv;

        // sum_i w_i w_i' over the sample, the only thing the beta block needs of
        // its regressors. Constant for the whole chain.
        arma::mat w_cross;

        // Variable selection. Only one of the two is ever engaged.
        std::optional<SsvsBlock> a_ssvs;
        std::optional<BvsBlock> a_bvs;
        arma::mat z_bvs;

        if (use_a)
        {
            diag_tt = arma::speye<arma::sp_mat>(tt, tt);
            diag_k = arma::eye(k, k);
            prior_a_mu = input.a_prior.mu;
            prior_a_vinv = input.a_prior.v_inv;
            a = input.initial.a;
            out.a = arma::mat(n_a, iterations);

            if (use_varsel)
            {
                out.a_lambda = arma::mat(n_a, iterations);

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

        if (use_beta)
        {
            diag_r = arma::eye(rank, rank);
            beta = input.initial.beta;
            w_cross = w_t * arma::trans(w_t);
            beta_mat = arma::reshape(beta, k_beta, rank);
            coint_v_inv = input.beta_prior.v_inv;
            coint_p_tau_inv = input.beta_prior.p_tau_inv;
            out.beta = arma::mat(n_beta, iterations);
        }

        // Error term
        const int post_u_sigma_df = input.u_sigma_prior.df + tt + rank;
        const arma::mat &prior_u_sigma_scale = input.u_sigma_prior.scale;
        arma::mat u_sigma_scale;
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

                // Block 1: Draw non-cointegration coefficients ----
                if (use_beta)
                {
                    // Update priors for alpha
                    // v^-1 (beta' P_tau^-1 beta) kron G^-1, the cointegration
                    // space prior of Koop, Leon-Gonzalez and Strachan (2010).
                    // The leading factor is the scalar shrinkage, not P_tau^-1
                    // again: as a matrix it conforms only when rank happens to
                    // equal k_beta, and the product of two symmetric matrices is
                    // not symmetric, so the posterior precision it built failed
                    // its Cholesky.
                    prior_a_vinv.submat(0, 0, n_alpha - 1, n_alpha - 1) = arma::kron(coint_v_inv * (arma::trans(beta_mat) * coint_p_tau_inv * beta_mat), u_sigma_inv);

                    // Update z matrix for alpha
                    if (a_bvs)
                    {
                        z_bvs.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta_mat) * w_t), diag_k);
                    }
                    else
                    {
                        z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta_mat) * w_t), diag_k);
                    }
                }

                if (a_bvs)
                {
                    z = z_bvs * a_bvs->lambda_diag;
                }

                // Update a
                //
                // The precision is symmetric, so z' D = (D z)' and one sparse
                // product serves both the posterior precision and its
                // right-hand side. Keeping the sparse operand on the left is
                // also what picks Armadillo's sparse-times-dense path rather
                // than promoting the whole block diagonal back to dense.
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
                    // trace(S U U') for the k x tt error matrix U. Contracting
                    // over k rather than forming the (k tt) quadratic form
                    // matters here and nowhere else: this closure runs twice per
                    // selected coefficient per draw, so it is the hottest
                    // arithmetic in the sampler.
                    bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta)
                              {
                    const arma::mat res = arma::reshape(y - z * theta, k, tt);
                    return -arma::accu((u_sigma_inv * res) % res) / 2; });
                }
            }

            // Block 2: Draw cointegration coefficients ----
            if (use_beta)
            {

                if (use_non_alpha)
                {
                    y_beta = y - z.cols(n_alpha, n_a - 1) * a.subvec(n_alpha, n_a - 1);
                }
                else
                {
                    y_beta = y;
                }

                // Reparameterise alpha
                alpha = arma::reshape(a.subvec(0, n_alpha - 1), k, rank);
                Alpha = alpha * arma::solve(arma::sqrtmat_sympd(arma::trans(alpha) * alpha), diag_r);

                // Update beta
                prior_beta_vinv = arma::kron(Alpha.t() * u_sigma_inv * Alpha, coint_v_inv * coint_p_tau_inv);
                post_beta_v = prior_beta_vinv + arma::kron(arma::trans(Alpha) * u_sigma_inv * Alpha, w_cross);
                Beta = draw_normal_precision(post_beta_v,
                                             arma::vectorise(w_t * arma::trans(arma::reshape(y_beta, k, tt)) *
                                                             u_sigma_inv * Alpha));
                Beta_mat = arma::reshape(Beta, k_beta, rank);

                // Final cointegration values. Only the product alpha beta' is
                // identified, so the draw is split between the two by the
                // normalisation below -- and both halves have to be carried
                // forward together.
                BB_sqrt = arma::sqrtmat_sympd(arma::trans(Beta_mat) * Beta_mat);
                beta_mat = Beta_mat * arma::solve(BB_sqrt, diag_r);

                alpha = Alpha * BB_sqrt;
                a.subvec(0, n_alpha - 1) = arma::vectorise(alpha);

                z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta_mat) * w_t), diag_k);
                u = arma::reshape(y - z * a, k, tt);
            }
            else
            {
                if (use_a)
                {
                    u = arma::reshape(y - z * a, k, tt);
                }
                else
                {
                    u = arma::reshape(y, k, tt);
                }
            }

            // Block 3: Draw inverse covariance matrix ----
            //
            // The cointegration space prior conditions alpha on Sigma, so it
            // contributes to Sigma's posterior in turn: the scale picks up
            // v^-1 alpha (beta' P_tau^-1 beta) alpha' and the degrees of freedom
            // pick up the r that goes with it. wishart() adds U U' itself.
            //
            // The prior scale is added. Left out, the inverse Wishart
            // prior has no influence on Sigma whatsoever and /priors/u_sigma/scale
            // is dead data that a caller has every reason to think is being used;
            // it is also inconsistent with counting r in the degrees of freedom,
            // which that branch does not do either.
            // Guarded, because a VEC of rank zero has no alpha to contribute.
            u_sigma_scale = prior_u_sigma_scale;
            if (use_beta)
            {
                u_sigma_scale += coint_v_inv * alpha *
                                 (beta_mat.t() * coint_p_tau_inv * beta_mat) * alpha.t();
            }
            u_sigma_inv = wishart(u, u_sigma_scale, post_u_sigma_df);

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

                if (use_beta)
                {
                    out.beta.col(draw_pos) = arma::vectorise(beta_mat);
                }

                out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
            }
        }

        reporter.finish();
        return out;
    }

    ForecastDraws VecNormalWishartSampler::forecast(const VecNormalWishartInput &input,
                                                    const VecNormalWishartDraws &coefficients,
                                                    Reporter &reporter) const
    {
        // A VEC and its level VAR are the same model in two parameterisations, so
        // the forecast is the VAR's: rewrite each draw in the level
        // parameterisation and hand the result to the sampler that already
        // simulates a level path. Simulating differences here instead would mean
        // a second copy of that recursion, and the VEC's own regressors cannot
        // drive it anyway -- their leading columns are loadings on the
        // cointegration space, not lags of the series a path is built from.
        //
        // `input.forecast.z` is therefore expected in the level VAR's layout,
        // the one vec_to_var_spec() describes: p + 1 lags of the endogenous
        // variables, one further block of unmodelled ones, and the deterministic
        // terms restricted to the cointegration space last. That is a real demand
        // on the caller, but the level history a forecast starts from cannot be
        // recovered from differenced regressors, so it has to arrive from
        // wherever the rest of the data is assembled.
        //
        // Nothing is validated here that the two calls below do not already
        // check: the conversion rejects draws that do not match the spec, and
        // the VAR's forecast rejects a `z` that does not match the converted
        // coefficients.
        VarNormalWishartInput var_input;
        var_input.spec = vec_to_var_spec(input.spec);
        var_input.forecast = input.forecast;

        const VarNormalWishartDraws var_coefficients =
            vec_to_var_coefficients(input.spec, coefficients);

        return VarNormalWishartSampler{}.forecast(var_input, var_coefficients, reporter);
    }

    arma::mat VecNormalWishartSampler::log_likelihood(const VecNormalWishartInput &input,
                                                      const VecNormalWishartDraws &coefficients) const
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
        arma::mat z = input.train.z;
        const arma::mat w_t = arma::trans(input.train.w);
        const int k_beta = static_cast<int>(input.spec.k_beta);
        const int rank = static_cast<int>(input.spec.rank);
        const bool use_a = z.n_cols > 0;
        const bool use_beta = rank > 0;
        const int n_alpha = static_cast<int>(input.spec.n_alpha());
        const arma::mat diag_k = arma::eye<arma::mat>(k, k);

        if (use_a && !coefficients.has_a())
        {
            throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
        }

        if (use_beta && !coefficients.has_beta())
        {
            throw std::invalid_argument("the model has beta coefficients but posterior draws of beta are missing");
        }

        const arma::uword draws = coefficients.iterations();
        const int tt = static_cast<int>(y.n_elem) / k;
        arma::mat loglik = arma::mat(draws, tt);

        // Calculate errors
        arma::mat u = arma::repmat(y, 1, draws);
        if (use_a)
        {
            if (use_beta)
            {
                for (arma::uword draw = 0; draw < draws; draw++)
                {
                    z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(arma::reshape(coefficients.beta.col(draw), k_beta, rank)) * w_t), diag_k);
                    u.col(draw) = y - z * coefficients.a.col(draw);
                }
            }
            else
            {
                u = u - z * coefficients.a;
            }
        }

        // Calculate log likelihood
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
