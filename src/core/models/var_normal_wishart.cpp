// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_wishart.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/constrained_var_path.h"
#include "core/algorithms/ssvs.h"
#include "core/algorithms/wishart.h"
#include "core/models/completion_support.h"
#include "core/models/constraint_support.h"
#include "core/models/forecast_states.h"
#include "core/models/model_support.h"
#include "core/models/shrinkage_support.h"
#include "core/models/steady_state_support.h"
#include "core/models/student_t_support.h"
#include "core/models/volatility_break_support.h"
#include "core/models/predictive_score.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace bayests
{

using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::covariance_root;
using core::draw_normal_precision;
using core::iid_block;
using core::IidBlock;
using core::SsvsBlock;
using core::ssvs_sweep;
using core::split_structural_coefficients;
using core::stacked_response;
using core::structural_inverse;
using core::require_forecast_regressors;
using core::update_forecast_lags;
using core::report_flat_selection_prior;

VarNormalWishartDraws VarNormalWishartSampler::draw_coefficients(const VarNormalWishartInput &input,
                                                                 Reporter &reporter) const
{
    input.validate();

    const int k = input.spec.k;
    const int iterations = input.spec.iterations;
    const int draws = input.spec.draws();

    // Not const: where the panel is not observed whole, every sweep completes it.
    arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;

    const int nparams = static_cast<int>(z.n_cols);
    const IidBlock iid = iid_block(input.spec, static_cast<arma::uword>(nparams));
    const bool use_a = nparams > 0;
    const int tt = static_cast<int>(y.n_elem) / k;

    const bool use_ssvs = input.spec.varsel == VarSelection::ssvs;
    const bool use_bvs = input.spec.varsel == VarSelection::bvs;
    const bool use_varsel = use_ssvs || use_bvs;

    VarNormalWishartDraws out;

    // Coefficients
    arma::vec a, prior_a_rhs;
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
        // The i.i.d. block, if the model has one: the restricted equations'
        // columns leave the regressors, their rows and columns leave the prior,
        // and the zeros go back when a draw is stored. Every one of these is the
        // identity when nothing is restricted.
        //
        // The prior's precision-weighted mean is formed once here rather than
        // once per draw. It was a constant inside the loop before, so the draws
        // are unchanged; reduced, it is the free part of V mu, which is what the
        // prior conditional on the restricted coefficients being zero asks for.
        prior_a_rhs = iid.elements(input.a_prior.v_inv * input.a_prior.mu);
        prior_a_vinv = iid.block(input.a_prior.v_inv);
        a = iid.elements(input.initial.a);
        z = iid.columns(z);
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
                report_flat_selection_prior(reporter, input.spec.varsel, "a", input.varsel_prior,
                                            input.a_prior.v_inv);
            }
        }
    }

    // Error term
    const int post_u_sigma_df = input.u_sigma_prior.df + tt;
    const arma::mat &prior_u_sigma_scale = input.u_sigma_prior.scale;
    arma::mat u_sigma_inv = input.initial.u_sigma_inv;
    out.u_sigma_inv = arma::mat(k * k, iterations);
    arma::mat u = arma::reshape(y, k, tt);

    // A panel not observed whole. Each sweep starts by drawing what was not
    // observed given the current coefficients and precision, then rebuilds the
    // lag block of the regressors from the completed panel, and everything
    // after that is the sweep of a complete panel, unchanged. A file without
    // constraints never enters any of this, so its draws are untouched.
    // The error precision of each soft group is drawn there too, after every
    // completion, from what the completed panel leaves its rows to explain.
    core::PanelCompletion completion(input.spec, input.train, use_a, input.constraints_prior,
                                     input.initial.constraints_inv);
    const bool complete = completion.active();
    completion.allocate(out.y, out.constraints_inv, static_cast<arma::uword>(iterations),
                        static_cast<arma::uword>(tt));

    // Start simulation
    // The adaptive prior on the coefficients, and how many sweeps found no
    // stationary draw and kept the one before; see shrinkage_support.h. Both
    // idle -- no random number, no rescaled prior -- unless /model asks.
    core::CoefficientShrinkage shrinkage(input.spec.shrinkage, input.a_shrinkage_prior, input.a_prior,
                                         input.initial.a_shrinkage, input.initial.a_local,
                                         input.initial.a_theta, input.spec.burnin);
    if (shrinkage.active())
    {
        out.a_shrinkage = arma::mat(shrinkage.groups(), iterations);
        if (core::shrinkage_has_local(input.spec.shrinkage))
        {
            out.a_local = arma::mat(nparams, iterations);
        }
        if (shrinkage.draws_theta())
        {
            out.a_theta = arma::mat(shrinkage.groups(), iterations);
        }
    }
    int unstationary_sweeps = 0;

    // The unconditional mean under the steady-state prior.
    arma::vec mu = input.initial.mu;
    if (input.spec.steady_state)
    {
        out.mu = arma::mat(static_cast<arma::uword>(k), iterations);
    }

    // A volatility break: every sweep runs on the rows from the break on divided
    // by s_t -- the VAR of the rescaled data, which given the scales is the whole
    // model -- and ends by drawing the scales from the residuals on the original
    // scale. Inactive, it holds ones, draws nothing and rescales nothing.
    core::VolatilityBreak volatility(input.spec, input.volatility_prior, input.initial.volatility,
                                     static_cast<arma::uword>(tt), static_cast<arma::uword>(k));
    // Student t errors: a mixing scale lambda_t per period on top, drawn after
    // the Wishart like the break, the two multiplying. Inactive, all ones.
    core::StudentTErrors tails(input.spec, input.student_t_prior, input.initial.student_t_df,
                               static_cast<arma::uword>(tt), static_cast<arma::uword>(k));
    const bool rescaled = volatility.active() || tails.active();
    const arma::vec y_original = rescaled ? y : arma::vec();
    const arma::mat z_original = rescaled ? z : arma::mat();
    if (volatility.active())
    {
        out.volatility = arma::mat(volatility.theta().n_elem, iterations);
    }
    if (tails.active())
    {
        out.student_t_df = arma::mat(1, iterations);
        out.student_t_lambda = arma::mat(static_cast<arma::uword>(tt), iterations);
    }
    // What every row is divided by, s_t sqrt(lambda_t), and the first period
    // that is: the break's, unless every period has a lambda of its own.
    const auto row_scale = [&]() {
        return tails.active() ? arma::vec(volatility.scales() % arma::sqrt(tails.lambda()))
                              : volatility.scales();
    };
    const arma::uword first_scaled = tails.active() ? 0 : volatility.start();

    for (int draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(draw + 1, draws);

        if (rescaled)
        {
            core::rescale_rows(y_original, z_original, row_scale(), static_cast<arma::uword>(k),
                               first_scaled, y, z);
            if (a_bvs)
            {
                z_bvs = z;
            }
            if (!use_a)
            {
                u = arma::reshape(y, k, tt);
            }
        }

        if (complete)
        {
            const arma::mat path =
                completion.complete(iid.scatter(a), core::covariance_of(u_sigma_inv, "VarNormalWishart"));
            y = arma::vectorise(path);
            if (use_a)
            {
                z = iid.columns(completion.regressors());
                if (a_bvs)
                {
                    z_bvs = z;
                }
            }
            else
            {
                u = path;
            }
        }

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
            if (input.spec.steady_state)
            {
                // The prior on the unconditional mean: the lags given mu, then mu
                // given the lags, the intercept they imply written into `a`.
                if (!core::draw_steady_state(a, mu, y, z, u_sigma_inv_diag, input.a_prior,
                                             input.mu_prior, static_cast<arma::uword>(k),
                                             static_cast<arma::uword>(input.spec.p),
                                             input.spec.stationary))
                {
                    unstationary_sweeps++;
                }
            }
            else
            {
                // An adaptive prior: the precision this draw is made under, rescaled
                // by the scales of the last sweep.
                if (shrinkage.active())
                {
                    prior_a_vinv = shrinkage.precision();
                    prior_a_rhs = shrinkage.rhs();
                }
                post_a_v = prior_a_vinv + arma::trans(dz) * z;
                if (input.spec.stationary)
                {
                    const arma::vec a_rhs = prior_a_rhs + arma::trans(dz) * y;
                    if (!core::draw_stationary(
                            a, [&]() { return draw_normal_precision(post_a_v, a_rhs); },
                            static_cast<arma::uword>(k), static_cast<arma::uword>(input.spec.p)))
                    {
                        unstationary_sweeps++;
                    }
                }
                else
                {
                    a = draw_normal_precision(post_a_v,
                                              prior_a_rhs + arma::trans(dz) * y);
                }
                shrinkage.update(a);
            }

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

        // Update the mixing scales and the volatility break, from the residuals
        // on the original scale: lambda given the residuals with the break
        // divided out, then the break given them with the new lambda divided out.
        if (rescaled)
        {
            arma::mat u_original = u;
            u_original.each_row() %= arma::trans(row_scale());
            if (tails.active())
            {
                arma::mat e = u_original;
                e.each_row() /= arma::trans(volatility.scales());
                tails.update(e, u_sigma_inv, draw < input.spec.burnin);
            }
            if (volatility.active())
            {
                if (tails.active())
                {
                    u_original.each_row() /= arma::trans(arma::sqrt(tails.lambda()));
                }
                volatility.update(u_original, u_sigma_inv, draw < input.spec.burnin);
            }
        }

        // Store draws
        if (input.spec.keeps(draw))
        {
            const int draw_pos = input.spec.kept_index(draw);
            if (input.spec.steady_state)
            {
                out.mu.col(draw_pos) = mu;
            }
            if (shrinkage.active())
            {
                out.a_shrinkage.col(draw_pos) = shrinkage.scale();
                if (core::shrinkage_has_local(input.spec.shrinkage))
                {
                    out.a_local.col(draw_pos) = shrinkage.local();
                }
                if (shrinkage.draws_theta())
                {
                    out.a_theta.col(draw_pos) = shrinkage.theta();
                }
            }
            if (use_a)
            {
                out.a.col(draw_pos) = iid.scatter(a);
                if (use_varsel)
                {
                    out.a_lambda.col(draw_pos) = a_ssvs ? a_ssvs->lambda : a_bvs->lambda;
                }
            }
            out.u_sigma_inv.col(draw_pos) = arma::vectorise(u_sigma_inv);
            if (volatility.active())
            {
                out.volatility.col(draw_pos) = volatility.theta();
            }
            if (tails.active())
            {
                out.student_t_df(0, draw_pos) = tails.df();
                out.student_t_lambda.col(draw_pos) = tails.lambda();
            }
            completion.store(out.y, out.constraints_inv, static_cast<arma::uword>(draw_pos), y);
        }
    }

    reporter.finish();
    volatility.report(reporter, input.spec.burnin);
    tails.report(reporter, input.spec.burnin);
    if (unstationary_sweeps > 0)
    {
        reporter.message(std::to_string(unstationary_sweeps) + " of " + std::to_string(draws) +
                         " sweeps found no stationary draw of the coefficients in " +
                         std::to_string(core::kStationaryTries) +
                         " tries and kept the one before: the posterior may put much of its mass "
                         "on explosive coefficients");
    }
    return out;
}

ForecastDraws VarNormalWishartSampler::forecast(const VarNormalWishartInput &input,
                                                const VarNormalWishartDraws &coefficients,
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
        throw std::invalid_argument("forecast regressors were supplied but posterior draws of a "
                                    "are missing");
    }

    // Counted off the posterior, not off x: the structural coefficients are the
    // last n_structural rows of a and have no column in x, so k * x.n_cols is
    // short by exactly that many.
    //
    // This is the path every VEC forecast reaches, converted to its level
    // parameterisation -- so leaving the split out here left a structural VEC
    // unable to forecast at all.
    const int nparams = (have_x || structural) ? static_cast<int>(coefficients.a.n_rows) : 0;
    const bool use_a = nparams > 0 && nparams > n_structural;

    arma::mat a = coefficients.a;
    const arma::mat a0 = split_structural_coefficients(input.spec, a, nparams);

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

    arma::mat fcst = arma::zeros<arma::mat>(h * k, draws);

    // A panel not observed whole: the lags a forecast starts from are the
    // draw's own completed panel, not what the host wrote into x, which for an
    // entry nothing observed is a placeholder.
    const bool complete = !input.train.constraints.empty();
    const arma::uword tt =
        complete ? stacked_response(input.train).n_elem / static_cast<arma::uword>(k) : 0;
    core::require_completed_panels(input, coefficients, use_a);
    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::mat error_root;

    // A scenario: the horizon is drawn whole, as the path the draw's VAR puts
    // on it conditioned on /data/forecast/constraints -- the completion step,
    // over h periods that start where the draw's panel ends. Only where a
    // scenario is given: an unconditional forecast keeps the recursion below,
    // and with it every draw it made before.
    const bool conditioned = !input.forecast.constraints.empty();

    // A volatility break: each horizon's error is drawn with the scale the
    // draw's path puts on it, periods tt .. tt + h - 1 counted from the start of
    // the sample, whether or not a scenario conditions the path.
    const bool scaled = input.spec.uses_volatility_break();
    const arma::uword tt_sample =
        scaled ? stacked_response(input.train).n_elem / static_cast<arma::uword>(k) : 0;
    if (scaled && coefficients.volatility.n_cols != draws)
    {
        throw std::invalid_argument("the model has a volatility break, and posterior draws of it "
                                    "(/posterior/volatility/coeffs) are missing or not one column "
                                    "per draw");
    }

    // Student t errors: each horizon's shock is drawn with a mixing scale of
    // its own from IG(nu / 2, nu / 2), nu the draw's.
    const bool tailed = input.spec.student_t;
    if (tailed && coefficients.student_t_df.n_cols != draws)
    {
        throw std::invalid_argument("the model has Student t errors, and posterior draws of their "
                                    "degrees of freedom (/posterior/student_t/df) are missing or "
                                    "not one column per draw");
    }
    if (tailed && conditioned)
    {
        throw std::invalid_argument("VarNormalWishart does not condition a forecast with Student t "
                                    "errors on /data/forecast/constraints yet");
    }

    // Calculate forecasts
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(draws));

        const arma::vec horizon_scale =
            scaled ? core::volatility_path(input.spec, coefficients.volatility.col(draw), tt_sample,
                                           static_cast<arma::uword>(h))
                   : arma::vec();

        if (conditioned)
        {
            fcst.col(draw) = arma::vectorise(core::conditioned_forecast(
                input, coefficients, draw, use_a, "VarNormalWishart", horizon_scale));
            continue;
        }

        // Once per draw: nothing in either depends on the horizon.
        const arma::mat a0_inv =
            structural ? structural_inverse(a0, draw, diag_k) : arma::mat();
        // The draw's coefficients as the k x n_x matrix they are. The SUR
        // spelling this replaced had them as a vector and paid for the
        // reshape implicitly, once per horizon, by widening z instead.
        const arma::mat a_draw =
            use_a ? arma::reshape(a.col(draw), k, x.n_cols) : arma::mat();

        // The error covariance factorised once per draw rather than once per
        // horizon: the precision is the same at every horizon, and the
        // factorisation draws nothing, so where it sits does not move a draw.
        error_root = covariance_root(arma::reshape(coefficients.u_sigma_inv.col(draw), k, k));

        // Every lag that reaches back into the sample, from this draw's panel:
        // horizon i's lag j for j > i. The ones inside the horizon are the
        // forecast's own and update_forecast_lags() writes them below.
        if (complete && use_a && p_larger_than_0)
        {
            core::start_from_panel(x, arma::reshape(coefficients.y.col(draw), k, tt),
                                   static_cast<arma::uword>(p));
        }

        for (int i = 0; i < h; i++)
        {
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
            if (scaled || tailed)
            {
                double shock_scale = scaled ? horizon_scale(static_cast<arma::uword>(i)) : 1.0;
                if (tailed)
                {
                    const double nu = coefficients.student_t_df(0, draw);
                    shock_scale *= std::sqrt(1.0 / arma::randg<double>(arma::distr_param(0.5 * nu, 2.0 / nu)));
                }
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) =
                    fcst.submat(i * k, draw, (i + 1) * k - 1, draw) +
                    shock_scale * (error_root * arma::randn(k));
            }
            else
            {
                fcst.submat(i * k, draw, (i + 1) * k - 1, draw) =
                    fcst.submat(i * k, draw, (i + 1) * k - 1, draw) + error_root * arma::randn(k);
            }

            // A_0 y_t = A_1 y_{t-1} + ... + u_t, so the inverse applies to the
            // whole right-hand side, signal and error alike.
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

namespace
{

/// The pointwise log likelihood of a complete panel, draws x tt: `y` stacked by
/// period, `z` its SUR regressors, and `first` the position of its first period
/// counted from the start of the training sample -- zero for the sample itself,
/// tt for a scored horizon -- which is what places a volatility break on it.
arma::mat gaussian_log_likelihood(const VarSpec &spec, const arma::vec &y, const arma::mat &z,
                                  const VarNormalWishartDraws &coefficients, const arma::uword first)
{
    const int k = spec.k;
    const bool use_a = z.n_cols > 0;
    const arma::uword draws = coefficients.iterations();
    const int tt = static_cast<int>(y.n_elem) / k;
    arma::mat loglik = arma::mat(draws, tt);

    const bool scaled = spec.uses_volatility_break();
    if (scaled && coefficients.volatility.n_cols != draws)
    {
        throw std::invalid_argument("the model has a volatility break, and posterior draws of it "
                                    "(/posterior/volatility/coeffs) are missing or not one column "
                                    "per draw");
    }
    const bool tailed = spec.student_t;
    if (tailed && coefficients.student_t_df.n_cols != draws)
    {
        throw std::invalid_argument("the model has Student t errors, and posterior draws of their "
                                    "degrees of freedom (/posterior/student_t/df) are missing or "
                                    "not one column per draw");
    }

    // Calculate errors
    arma::mat u = arma::repmat(y, 1, draws);
    if (use_a)
    {
        u = u - z * coefficients.a;
    }

    // Calculate log likelihood. Under a volatility break period t's covariance
    // is s_t^2 Sigma: the log determinant gains -k log s_t, the quadratic form
    // is divided by s_t^2.
    const double part_a = -k * std::log(2 * arma::datum::pi) / 2;
    arma::mat u_sigma_inv;
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        u_sigma_inv = arma::reshape(coefficients.u_sigma_inv.col(draw), k, k);
        const double part_b = core::half_log_det_precision(u_sigma_inv);
        const arma::vec s =
            scaled ? core::volatility_path(spec, coefficients.volatility.col(draw), first,
                                           static_cast<arma::uword>(tt))
                   : arma::vec();
        for (int i = 0; i < tt; i++)
        {
            const double part_c = -arma::as_scalar(arma::trans(u.submat(i * k, draw, (i + 1) * k - 1, draw)) * u_sigma_inv * u.submat(i * k, draw, (i + 1) * k - 1, draw)) / 2;
            if (tailed)
            {
                // The multivariate t with scale s_t^2 Sigma: lambda_t integrated
                // out, so the density of the data rather than of one draw's mix.
                const double si = scaled ? s(static_cast<arma::uword>(i)) : 1.0;
                loglik(draw, i) = core::student_t_log_density(
                    -2.0 * part_c / (si * si), -2.0 * part_b + 2.0 * k * std::log(si),
                    static_cast<double>(k), coefficients.student_t_df(0, draw));
            }
            else if (scaled)
            {
                const double si = s(static_cast<arma::uword>(i));
                loglik(draw, i) = part_a + part_b - k * std::log(si) + part_c / (si * si);
            }
            else
            {
                loglik(draw, i) = part_a + part_b + part_c;
            }
        }
    }

    return loglik;
}

/// The likelihood and the score run from draws without validate() in between,
/// so the one combination they would get wrong is refused here as well: the
/// completion step's density has one covariance for every period.
void require_no_break_with_constraints(const VarNormalWishartInput &input)
{
    if (input.spec.uses_volatility_break() || input.spec.student_t)
    {
        throw std::invalid_argument(
            "VarNormalWishart does not combine /model/volatility_break or /model/student_t with "
            "/data/train/constraints or /data/test/constraints");
    }
}

} // namespace

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

    // A panel not observed whole: the density of what was observed, with what
    // was not integrated out -- never the density of one draw's completion,
    // which is not a likelihood of the data. Period t is the density of the
    // rows ending in t given those ending before it.
    if (!input.train.constraints.empty())
    {
        require_no_break_with_constraints(input);
        return core::observed_log_likelihood(input, coefficients, "VarNormalWishart");
    }

    return gaussian_log_likelihood(input.spec, y, z, coefficients, 0);
}

arma::mat VarNormalWishartSampler::predictive_log_density(const VarNormalWishartInput &input,
                                                        const VarNormalWishartDraws &coefficients) const
{
    core::require_scorable(input.spec, "VarNormalWishart");
    const arma::uword periods = core::scored_horizons(input.test.y, input.spec);
    const arma::mat x =
        core::realised_regressors(input.forecast.x, input.test.y, input.spec.k, input.spec.p);

    // A sample not observed whole, or a horizon realised in part: the density
    // of what the horizon realised under the path each draw puts on it, from
    // where that draw's panel ends, period by period -- the completion step's
    // log density. For a horizon realised whole it is the same pointwise log
    // likelihood as below, computed another way; the recursion below stays for
    // every file without constraints, so its scores do not move.
    if (!input.train.constraints.empty() || !input.test.constraints.empty())
    {
        require_no_break_with_constraints(input);
        return core::score_from_panel(input, coefficients, periods, "VarNormalWishart");
    }

    // The scored periods as a sample of their own. Nothing in this model moves
    // over the horizon, so a draw describes period T + i exactly as it
    // describes the sample, and the density of the realised values under it is
    // the pointwise log likelihood over that sample -- the same expression, the
    // same code, a different tt.
    //
    // Only the spec and the three members below are filled, which is everything
    // log_likelihood() reads. unit.predictive_score pins the two against each
    // other on a sample this is handed back verbatim, so a density that grew a
    // dependency would fail there rather than quietly score against a default.
    VarNormalWishartInput scored;
    scored.spec = input.spec;
    scored.train.y = input.test.y.head_rows(periods);
    scored.train.x = x;
    scored.train.z = core::sur_regressors(x, input.spec.k);

    // Under a volatility break the horizon's periods are tt and on, which is
    // where the path puts their scales; without one this is log_likelihood().
    if (input.spec.uses_volatility_break())
    {
        if (coefficients.u_sigma_inv.n_elem == 0)
        {
            throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
        }
        const arma::uword tt_sample =
            stacked_response(input.train).n_elem / static_cast<arma::uword>(input.spec.k);
        return gaussian_log_likelihood(scored.spec, stacked_response(scored.train), scored.train.z,
                                       coefficients, tt_sample);
    }
    return log_likelihood(scored, coefficients);
}

} // namespace bayests
