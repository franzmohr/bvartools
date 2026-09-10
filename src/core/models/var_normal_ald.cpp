// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_ald.h"

#include "core/algorithms/bvs.h"
#include "core/models/ald_support.h"
#include "core/models/model_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::AldShape;
using core::ald_log_density;
using core::ald_shape;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_ald_scale;
using core::draw_ald_weights;
using core::draw_normal_precision;
using core::stacked_response;

VarNormalAldDraws VarNormalAldSampler::draw_coefficients(const VarNormalAldInput &input,
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

    // The two constants the quantile enters through. At q = 0.5 the skew is
    // zero and tau2 is eight, which is the symmetric case the whole model
    // collapses to.
    const AldShape shape = ald_shape(input.spec.quantile);

    // Only BVS is implemented for this model.
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_bvs;

    VarNormalAldDraws out;

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

    // The error term, and the latent scales it is a mixture over. `w` is
    // tt x k, one column per equation, so vectorise(trans(w)) is period-major
    // and lines up with the stacked response.
    arma::mat u = arma::reshape(y, k, tt);
    arma::mat w = input.initial.w;
    arma::vec u_scale = input.initial.u_scale;

    // The response the coefficient block regresses on: the data less the skew
    // the mixture carries. At the median this subtracts nothing.
    arma::vec y_adjusted = y - shape.theta * arma::vectorise(arma::trans(w));

    // Sigma is diagonal in every period -- there is no covariance block for this
    // model -- so the sparsity structure never changes and only the diagonal is
    // ever written after this point.
    arma::sp_mat u_sigma_inv_diag = arma::eye<arma::sp_mat>(k * tt, k * tt);
    u_sigma_inv_diag.diag() =
        1 / (shape.tau2 * arma::vectorise(arma::trans(w)) %
             arma::repmat(u_scale, tt, 1));

    const arma::vec u_scale_post_shape =
        input.u_scale_prior.shape + static_cast<double>(tt) * 3.0 / 2.0;
    const arma::vec &u_scale_prior_rate = input.u_scale_prior.rate;

    out.u_scale = arma::mat(k, iterations);
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
            a = draw_normal_precision(a_post_v, a_prior_vinv * a_prior_mu +
                                                    arma::trans(z) * u_sigma_inv_diag * y_adjusted);

            if (use_bvs)
            {
                // The skew has to be inside the residual here as well. Scoring
                // an inclusion against y - z*theta rather than
                // y - theta*w - z*theta would select for the median whatever
                // quantile was asked for, and at q = 0.5 nothing would show.
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::vec res = y_adjusted - z_bvs * a_bvs->lambda_diag * theta;
                    return -arma::as_scalar(arma::trans(res) * u_sigma_inv_diag * res) / 2;
                });
                z = z_bvs * a_bvs->lambda_diag;
            }

            u = arma::reshape(y - z * a, k, tt);
        }
        else
        {
            u = arma::reshape(y, k, tt);
        }

        // Update the latent scales, one per observation. Their full conditional
        // is generalised inverse Gaussian at index 1/2.
        w = draw_ald_weights(u, u_scale, shape);

        // Update the scale of the asymmetric Laplace, one per equation.
        u_scale = draw_ald_scale(u, w, shape, u_scale_post_shape, u_scale_prior_rate);

        // Rebuild what the next sweep regresses against.
        const arma::vec w_stacked = arma::vectorise(arma::trans(w));
        y_adjusted = y - shape.theta * w_stacked;
        u_sigma_inv_diag.diag() =
            1 / (shape.tau2 * w_stacked % arma::repmat(u_scale, tt, 1));

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

            out.u_scale.col(draw_pos) = u_scale;
            out.u_omega_inv.col(draw_pos) = u_sigma_inv_diag.diag();

            for (int i = 0; i < tt; i++)
            {
                out.u_sigma_inv.submat(i * k * k, draw_pos, (i + 1) * k * k - 1, draw_pos) =
                    arma::vectorise(arma::mat(u_sigma_inv_diag.submat(
                        i * k, i * k, (i + 1) * k - 1, (i + 1) * k - 1)));
            }
        }
    }

    reporter.finish();
    return out;
}

ForecastDraws VarNormalAldSampler::forecast(const VarNormalAldInput &input,
                                            const VarNormalAldDraws &draws,
                                            Reporter &reporter) const
{
    (void)input;
    (void)draws;
    (void)reporter;

    throw std::invalid_argument(
        "a quantile regression model does not forecast: the h step ahead quantile is not the "
        "quantile of the iterated one step ahead quantiles, so iterating this model forward "
        "would produce a path that cannot be read as a quantile of anything");
}

arma::mat VarNormalAldSampler::log_likelihood(const VarNormalAldInput &input,
                                              const VarNormalAldDraws &coefficients) const
{
    const int k = input.spec.k;

    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (coefficients.u_scale.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of the asymmetric Laplace scale are missing");
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
    const double q = input.spec.quantile;

    // Errors, one column per draw.
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        u = u - z * coefficients.a;
    }

    // The asymmetric Laplace density itself, which is closed form and marginal
    // of the latent scales -- so no state is conditioned on and no determinant
    // is recomputed. A period's contribution is the sum over the k equations,
    // each under its own scale.
    arma::mat loglik(draws, tt);
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        const arma::vec scale = coefficients.u_scale.col(draw);
        for (int i = 0; i < tt; i++)
        {
            double total = 0.0;
            for (int j = 0; j < k; j++)
            {
                total += ald_log_density(u(static_cast<arma::uword>(i * k + j), draw), scale(j), q);
            }
            loglik(draw, static_cast<arma::uword>(i)) = total;
        }
    }

    return loglik;
}

} // namespace bayests
