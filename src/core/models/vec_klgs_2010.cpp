// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/vec_klgs_2010.h"

#include "bayests/var_normal_wishart.h"
#include "bayests/vec_to_var.h"
#include "core/algorithms/wishart.h"
#include "core/models/model_support.h"
#include "core/models/vec_support.h"

#include <cmath>
#include <stdexcept>

namespace bayests
{

namespace
{

using core::draw_normal_precision;
using core::reparameterise_alpha;
using core::stacked_response;

/// The response as the compact form wants it: k x tt, one period per column.
///
/// Built through stacked_response() rather than by transposing `train.y`,
/// because `y` is allowed to arrive already stacked -- a single row or a single
/// column of vec(y') is how the HDF5 files store it -- and only the stacked
/// vector is the same object in all three cases.
arma::mat response_by_period(const TrainData &train, const int k, const int tt)
{
    return arma::reshape(stacked_response(train), k, tt);
}

/// The transposed design matrix, n_design x tt, with the compact regressors
/// already in place below the error correction block.
///
/// The leading `rank` rows hold beta' w_t and are rewritten on every draw; the
/// rest is data and is written once. Kept transposed because everything that
/// reads it wants it that way -- the coefficient matrix multiplies it from the
/// left to give the k x tt errors, and its Gram product is the same either way.
arma::mat design_by_period(const arma::mat &x_t, const int n_design, const int rank, const int tt)
{
    arma::mat design_t(n_design, tt, arma::fill::zeros);
    if (n_design > rank)
    {
        design_t.rows(rank, n_design - 1) = x_t;
    }
    return design_t;
}

/// The same posterior in the layout the shared machinery is written against.
///
/// vec_to_var_coefficients() takes a VecNormalWishartDraws, and this model's
/// output is that one without the inclusion indicators -- same coefficient
/// order, same beta, same precision -- so the conversion is a move rather than a
/// transformation. Worth the copy: the alternative is a second overload of the
/// VEC-to-VAR change of basis that differs from the first in nothing.
VecNormalWishartDraws as_wishart_draws(const VecKlgs2010Draws &draws)
{
    VecNormalWishartDraws out;
    out.a = draws.a;
    out.beta = draws.beta;
    out.u_sigma_inv = draws.u_sigma_inv;
    return out;
}

} // namespace

VecKlgs2010Draws VecKlgs2010Sampler::draw_coefficients(const VecKlgs2010Input &input,
                                                       Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int iterations = input.spec.iterations;
    const int burnin = input.spec.burnin;
    const int draws = input.spec.draws();

    const int rank = input.spec.rank;
    const int k_beta = input.spec.k_beta;
    const int n_alpha = input.spec.n_alpha();
    const int n_beta = input.spec.n_beta();
    const int n_a = input.n_a();
    const int n_design = input.n_design();

    const bool use_a = n_a > 0;
    const bool use_beta = rank > 0;
    const bool use_non_alpha = n_design > rank;

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat y_t = response_by_period(input.train, k, tt);
    const arma::mat w_t = arma::trans(input.train.w);
    const arma::mat x_t = arma::trans(input.train.x);

    VecKlgs2010Draws out;

    // Coefficients (non-cointegration). `a_mat` is the k x n_design reading of
    // `a`: vec of it is `a` itself, so the two are the same numbers and neither
    // is derived from the other more than once per draw.
    arma::vec a, prior_a_mu;
    arma::mat prior_a_vinv, a_mat, post_a_v, design_t;

    // Coefficients (cointegration)
    arma::mat alpha, Alpha, beta_mat, Beta_mat, BB_sqrt, diag_r;
    arma::mat prior_beta_vinv, post_beta_v, alpha_sigma_alpha;
    // Initialised, not merely declared: block 3 reads it whenever the model
    // has a cointegration relation, and a rank-zero VEC never assigns it.
    double coint_v_inv = 0.0;
    arma::mat coint_p_tau_inv;

    // sum_t w_t w_t' over the sample, the only thing the beta block needs of
    // its regressors. Constant for the whole chain.
    arma::mat w_cross;

    if (use_a)
    {
        prior_a_mu = input.a_prior.mu;
        prior_a_vinv = input.a_prior.v_inv;
        a = input.initial.a;
        a_mat = arma::reshape(a, k, n_design);
        design_t = design_by_period(x_t, n_design, rank, tt);
        out.a = arma::mat(n_a, iterations);
    }

    if (use_beta)
    {
        diag_r = arma::eye(rank, rank);
        beta_mat = arma::reshape(input.initial.beta, k_beta, rank);
        w_cross = w_t * arma::trans(w_t);
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
    arma::mat u = y_t;

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
                // Update priors for alpha
                // v^-1 (beta' P_tau^-1 beta) kron G^-1, the cointegration
                // space prior of Koop, Leon-Gonzalez and Strachan (2010). The
                // leading factor is the scalar shrinkage, not P_tau^-1 again:
                // as a matrix it conforms only when rank happens to equal
                // k_beta, and the product of two symmetric matrices is not
                // symmetric, so the posterior precision it built would fail its
                // Cholesky.
                prior_a_vinv.submat(0, 0, n_alpha - 1, n_alpha - 1) =
                    arma::kron(coint_v_inv * (arma::trans(beta_mat) * coint_p_tau_inv * beta_mat),
                               u_sigma_inv);

                // Update the loadings' regressors, which are beta' w_t rather
                // than data. One product for the whole sample, r x tt.
                design_t.rows(0, rank - 1) = arma::trans(beta_mat) * w_t;
            }

            // Update a
            //
            // This is where the compact form pays. The SUR posterior precision
            // is z' kron(I_tt, Sigma^-1) z with z = kron(W_x, I_k), and that
            // Kronecker structure factors straight out of the product:
            //
            //     kron(W_x, I_k)' kron(I_tt, Sigma^-1) kron(W_x, I_k)
            //         = kron(W_x' W_x, Sigma^-1),
            //
            // by kron(A, B) kron(C, D) = kron(AC, BD) applied twice. The Gram
            // product left to form is n_design square over tt periods instead
            // of n_a square over tt k rows, and no (tt k) x n_a matrix is built
            // at all. The right-hand side collapses the same way, to
            // vec(Sigma^-1 Y' W_x).
            post_a_v = prior_a_vinv + arma::kron(design_t * arma::trans(design_t), u_sigma_inv);
            a = draw_normal_precision(
                post_a_v, prior_a_vinv * prior_a_mu +
                              arma::vectorise(u_sigma_inv * (y_t * arma::trans(design_t))));
            a_mat = arma::reshape(a, k, n_design);
        }

        // Block 2: Draw cointegration coefficients ----
        //
        // Unchanged from the SUR sampler, and not for want of trying: beta's
        // own regressors are kron(alpha, w_t'), which varies with t through w_t
        // and so has no structure to factor out. This block is compact in both.
        if (use_beta)
        {
            const arma::mat y_beta_t =
                use_non_alpha ? arma::mat(y_t - a_mat.cols(rank, n_design - 1) * x_t) : y_t;

            // Reparameterise alpha
            alpha = a_mat.cols(0, rank - 1);
            Alpha = reparameterise_alpha(alpha, diag_r);

            // Update beta
            alpha_sigma_alpha = arma::trans(Alpha) * u_sigma_inv * Alpha;
            prior_beta_vinv = arma::kron(alpha_sigma_alpha, coint_v_inv * coint_p_tau_inv);
            post_beta_v = prior_beta_vinv + arma::kron(alpha_sigma_alpha, w_cross);
            Beta_mat = arma::reshape(
                draw_normal_precision(post_beta_v, arma::vectorise(w_t * arma::trans(y_beta_t) *
                                                                   u_sigma_inv * Alpha)),
                k_beta, rank);

            // Final cointegration values. Only the product alpha beta' is
            // identified, so the draw is split between the two by the
            // normalisation below -- and both halves have to be carried
            // forward together.
            BB_sqrt = arma::sqrtmat_sympd(arma::trans(Beta_mat) * Beta_mat);
            beta_mat = Beta_mat * arma::solve(BB_sqrt, diag_r);

            alpha = Alpha * BB_sqrt;
            a_mat.cols(0, rank - 1) = alpha;
            a.subvec(0, n_alpha - 1) = arma::vectorise(alpha);

            design_t.rows(0, rank - 1) = arma::trans(beta_mat) * w_t;
        }

        u = use_a ? arma::mat(y_t - a_mat * design_t) : y_t;

        // Block 3: Draw inverse covariance matrix ----
        //
        // The cointegration space prior conditions alpha on Sigma, so it
        // contributes to Sigma's posterior in turn: the scale picks up
        // v^-1 alpha (beta' P_tau^-1 beta) alpha' and the degrees of freedom
        // pick up the r that goes with it. wishart() adds U U' itself.
        //
        // The prior scale is added, and the rank counted, exactly as
        // VecNormalWishart does and unlike bvartools' .simulation_klgs2010,
        // which overwrites the prior scale with the cointegration term and
        // leaves /priors/u_sigma/scale dead data that a caller has every reason
        // to think is being used.
        //
        // Guarded, because a VEC of rank zero has no alpha to contribute.
        u_sigma_scale = prior_u_sigma_scale;
        if (use_beta)
        {
            u_sigma_scale += coint_v_inv * alpha *
                             (arma::trans(beta_mat) * coint_p_tau_inv * beta_mat) *
                             arma::trans(alpha);
        }
        u_sigma_inv = wishart(u, u_sigma_scale, post_u_sigma_df);

        // Store draws
        if (draw >= burnin)
        {
            const int draw_pos = draw - burnin;

            if (use_a)
            {
                out.a.col(draw_pos) = a;
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

ForecastDraws VecKlgs2010Sampler::forecast(const VecKlgs2010Input &input,
                                           const VecKlgs2010Draws &coefficients,
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
    // `input.forecast.z` is therefore expected in the level VAR's layout, the
    // one vec_to_var_spec() describes, and stays a SUR matrix even though this
    // model's training regressors are compact: what reads it is
    // VarNormalWishartSampler, which is written against `z` and is shared with
    // every other VEC.
    //
    // Nothing is validated here that the two calls below do not already
    // check: the conversion rejects draws that do not match the spec, and
    // the VAR's forecast rejects a `z` that does not match the converted
    // coefficients.
    VarNormalWishartInput var_input;
    var_input.spec = vec_to_var_spec(input.spec);
    var_input.forecast = input.forecast;

    const VarNormalWishartDraws var_coefficients =
        vec_to_var_coefficients(input.spec, as_wishart_draws(coefficients));

    return VarNormalWishartSampler{}.forecast(var_input, var_coefficients, reporter);
}

arma::mat VecKlgs2010Sampler::log_likelihood(const VecKlgs2010Input &input,
                                             const VecKlgs2010Draws &coefficients) const
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

    const int rank = input.spec.rank;
    const int k_beta = input.spec.k_beta;
    const int n_design = input.n_design();
    const bool use_a = input.use_a();
    const bool use_beta = input.use_beta();

    if (use_a && !coefficients.has_a())
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }

    if (use_beta && !coefficients.has_beta())
    {
        throw std::invalid_argument("the model has beta coefficients but posterior draws of beta are missing");
    }

    const int tt = static_cast<int>(input.train.periods(k));
    const arma::mat y_t = response_by_period(input.train, k, tt);
    const arma::mat w_t = arma::trans(input.train.w);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);

    arma::mat design_t;
    if (use_a)
    {
        design_t = design_by_period(arma::trans(input.train.x), n_design, rank, tt);
    }

    const arma::uword draws = coefficients.iterations();
    arma::mat loglik(draws, tt);

    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u = y_t;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        // The loadings' regressors are a function of the draw, so the design
        // matrix is rebuilt from the beta that was stored beside it.
        if (use_beta)
        {
            design_t.rows(0, rank - 1) =
                arma::trans(arma::reshape(coefficients.beta.col(draw), k_beta, rank)) * w_t;
        }
        if (use_a)
        {
            u = y_t - arma::reshape(coefficients.a.col(draw), k, n_design) * design_t;
        }

        const arma::mat u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, k);
        const double part_b = -std::log(arma::det(arma::solve(u_sigma_inv, diag_k))) / 2;
        for (int i = 0; i < tt; i++)
        {
            const double part_c =
                -arma::as_scalar(arma::trans(u.col(i)) * u_sigma_inv * u.col(i)) / 2;
            loglik(draw, i) = part_a + part_b + part_c;
        }
    }

    return loglik;
}

} // namespace bayests
