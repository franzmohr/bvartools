// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_normal_gamma.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/bvs.h"
#include "core/algorithms/ssvs.h"
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
using core::reparameterise_alpha;
using core::ssvs_sweep;
using core::SsvsBlock;
using core::stacked_response;

VecNormalGammaDraws VecNormalGammaSampler::draw_coefficients(const VecNormalGammaInput &input,
                                                             Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
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

    const bool use_psi = input.use_psi();
    const int n_psi = input.spec.n_psi();

    const bool use_ssvs = input.spec.varsel == VarSelection::ssvs;
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_ssvs || use_bvs;

    VecNormalGammaDraws out;

    // Coefficients (non-cointegration)
    arma::vec a, a_prior_mu;
    arma::mat a_prior_vinv, a_post_v, dz, diag_k;

    // The error precision is the same in every period, so the block diagonal the
    // SUR form needs is kron(I_tt, u_sigma_inv): tt identical k x k blocks, and
    // hence a density of 1/tt. Same choice, for the same reason, as in
    // vec_normal_wishart.cpp.
    arma::sp_mat diag_tt, u_sigma_inv_diag;

    // Coefficients (cointegration)
    arma::mat alpha, Alpha, beta_mat, Beta_mat, BB_sqrt, diag_r;
    double coint_v_inv = 0.0;
    arma::mat coint_p_tau_inv, prior_beta_vinv, post_beta_v;
    arma::vec y_beta, Beta;

    // sum_i w_i w_i' over the sample, the only thing the beta block needs of its
    // regressors. Constant for the whole chain.
    arma::mat w_cross;

    // Variable selection. Only one of the two is ever engaged.
    std::optional<SsvsBlock> a_ssvs;
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;

    if (use_a)
    {
        diag_tt = arma::speye<arma::sp_mat>(tt, tt);
        diag_k = arma::eye(k, k);
        a_prior_mu = input.a_prior.mu;
        a_prior_vinv = input.a_prior.v_inv;
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
        w_cross = w_t * arma::trans(w_t);
        beta_mat = arma::reshape(input.initial.beta, k_beta, rank);
        coint_v_inv = input.beta_prior.v_inv;
        coint_p_tau_inv = input.beta_prior.p_tau_inv;
        out.beta = arma::mat(n_beta, iterations);
    }

    // Covariance block
    arma::vec psi, psi_prior_mu, psi_y;
    arma::mat Psi, Psi_lambda, psi_prior_vinv, psi_post_v, psi_z, dpsi_z, psi_u_omega_inv;
    arma::sp_mat psi_u_omega_inv_diag;

    std::optional<SsvsBlock> psi_ssvs;
    std::optional<BvsBlock> psi_bvs;
    arma::mat psi_z_bvs;

    if (use_psi)
    {
        psi_prior_mu = input.psi_prior.mu;
        psi_prior_vinv = input.psi_prior.v_inv;
        psi = input.initial.psi;
        psi_z = arma::zeros<arma::mat>(tt * (k - 1), n_psi);
        Psi = arma::eye<arma::mat>(k, k);
        Psi_lambda = arma::eye<arma::mat>(k, k);
        fill_strict_lower_triangle(Psi, psi);
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
                // zero matrix, exactly as in var_normal_gamma.cpp -- what the
                // BVS residuals are measured against is a modelling decision
                // rather than something to change while porting.
                psi_z_bvs = psi_z;
                psi_bvs.emplace(input.initial.psi_lambda, input.psi_varsel_prior);
            }
        }
    }

    // Error term
    const arma::vec u_sigma_post_shape = input.u_sigma_prior.shape + tt * 0.5;
    const arma::vec &u_sigma_prior_rate = input.u_sigma_prior.rate;
    // The chain starts from Psi' Omega^-1 Psi, with Psi unpacked above from the
    // k(k-1)/2 free elements /initial/psi carries -- the file stores the vector,
    // never the matrix. var_normal_gamma.cpp instead takes the stored precision
    // as Sigma^-1 directly, which leaves that vector unread: the psi block
    // overwrites it before anything looks at it, so the covariance block's
    // starting value has no effect there at all. Deriving it costs one product
    // before the loop and makes the datum mean what it says. It moves the first
    // draw of `a` and nothing after burn-in.
    arma::mat sse;
    arma::mat u_omega_inv = input.initial.u_sigma_inv;
    arma::mat u_sigma_inv = use_psi ? arma::mat(arma::trans(Psi) * u_omega_inv * Psi)
                                    : u_omega_inv;
    arma::mat u = arma::reshape(y, k, tt);

    out.u_omega_inv = arma::mat(k, iterations);
    out.u_sigma_inv = arma::mat(k * k, iterations);

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
                // v^-1 (beta' P_tau^-1 beta) kron G^-1, the cointegration space
                // prior of Koop, Leon-Gonzalez and Strachan (2010). The G it
                // conditions on is the error precision -- constant here, so
                // simply u_sigma_inv; VecNormalStochvol has to average over the
                // sample for it. The leading factor is the scalar shrinkage, not
                // P_tau^-1 again: see the same line in vec_normal_wishart.cpp.
                a_prior_vinv.submat(0, 0, n_alpha - 1, n_alpha - 1) =
                    arma::kron(coint_v_inv * (arma::trans(beta_mat) * coint_p_tau_inv * beta_mat),
                               u_sigma_inv);

                // Update z matrix for alpha
                fill_z_alpha_constant(a_bvs ? z_bvs : z, beta_mat, w_t, n_alpha, diag_k);
            }

            if (a_bvs)
            {
                z = z_bvs * a_bvs->lambda_diag;
            }

            // Update a. The precision is symmetric, so z' D = (D z)' and one
            // sparse product serves both the posterior precision and its
            // right-hand side.
            dz = u_sigma_inv_diag * z;
            a_post_v = a_prior_vinv + arma::trans(dz) * z;
            a = draw_normal_precision(a_post_v, a_prior_vinv * a_prior_mu + arma::trans(dz) * y);

            if (a_ssvs)
            {
                ssvs_sweep(*a_ssvs, a, a_prior_vinv);
            }

            if (a_bvs)
            {
                z = z_bvs;
                // res' kron(I_tt, S) res is sum_t u_t' S u_t, that is
                // trace(S U U') for the k x tt error matrix U.
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::mat res = arma::reshape(y - z * theta, k, tt);
                    return -arma::accu((u_sigma_inv * res) % res) / 2;
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

            // Update beta. The precision is the same in every period, so
            // sum_t z_t' S z_t collapses to kron(Alpha' S Alpha, sum_t w_t w_t')
            // and the whole sample enters through w_cross.
            prior_beta_vinv = arma::kron(Alpha.t() * u_sigma_inv * Alpha,
                                         coint_v_inv * coint_p_tau_inv);
            post_beta_v = prior_beta_vinv + arma::kron(arma::trans(Alpha) * u_sigma_inv * Alpha,
                                                       w_cross);
            Beta = draw_normal_precision(
                post_beta_v, arma::vectorise(w_t * arma::trans(arma::reshape(y_beta, k, tt)) *
                                             u_sigma_inv * Alpha));
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
                fill_strict_lower_triangle(Psi_lambda,
                                           psi_ssvs ? psi_ssvs->lambda : psi_bvs->lambda);
            }
            u = Psi * u;
        }

        // Block 4: Draw the error precisions ----
        //
        // No cointegration term is added, unlike VecNormalWishart, whose scale
        // picks up v^-1 alpha (beta' P_tau^-1 beta) alpha' and whose degrees of
        // freedom pick up the rank that goes with it. The prior conditions alpha
        // on the error precision either way, so a term is owed back; independent
        // gammas have no conjugate form for it, and bvartools' .bvecalg adds
        // none. Leaving it out is therefore the documented behaviour rather than
        // an oversight, but it does mean the two VECs treat the same prior
        // slightly differently.
        sse = u * u.t();
        for (int i = 0; i < k; i++)
        {
            u_omega_inv(i, i) = arma::randg<double>(arma::distr_param(
                u_sigma_post_shape(i),
                1 / arma::as_scalar(u_sigma_prior_rate(i) + sse(i, i) * 0.5)));
        }
        u_sigma_inv = use_psi ? arma::mat(arma::trans(Psi) * u_omega_inv * Psi) : u_omega_inv;

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

ForecastDraws VecNormalGammaSampler::forecast(const VecNormalGammaInput &input,
                                              const VecNormalGammaDraws &coefficients,
                                              Reporter &reporter) const
{
    // A VEC and its level VAR are the same model in two parameterisations, so
    // the forecast is the VAR's -- see VecNormalWishartSampler::forecast() for
    // why that is worth doing and what it demands of `input.forecast.z`.
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

arma::mat VecNormalGammaSampler::log_likelihood(const VecNormalGammaInput &input,
                                                const VecNormalGammaDraws &coefficients) const
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

    // Calculate errors. The loadings' regressors are a function of beta, so with
    // a cointegration relation they are rebuilt once per draw.
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

    // Calculate log likelihood
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, k);
        const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
        for (int i = 0; i < tt; i++)
        {
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
