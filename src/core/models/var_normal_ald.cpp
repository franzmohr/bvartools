// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/var_normal_ald.h"

#include "core/algorithms/bvs.h"
#include "core/algorithms/triangular_packing.h"
#include "core/models/ald_support.h"
#include "core/models/model_support.h"
#include "core/models/quantile_grid.h"

#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

namespace bayests
{

namespace
{

/// Progress over a whole quantile grid: each quantile's run reports its own
/// draws, and this maps them onto one count, so a host sees a single bar.
class GridReporter final : public Reporter
{
public:
    GridReporter(Reporter &parent, const long long offset, const long long total)
        : parent_(parent), offset_(offset), total_(total)
    {
    }
    void message(const std::string &text) override { parent_.message(text); }
    void progress(const long long done, const long long) override
    {
        parent_.progress(offset_ + done, total_);
    }
    void finish() override {}
    void check_interrupt() override { parent_.check_interrupt(); }

private:
    Reporter &parent_;
    long long offset_;
    long long total_;
};

/// A uniform strictly inside (0, 1), so that a tail of the grid is never asked
/// for the level zero.
double open_uniform()
{
    double u = 0.0;
    while (!(u > 0.0 && u < 1.0))
    {
        u = arma::randu<double>();
    }
    return u;
}

/// One draw's coefficients at every level of the grid: the k x n_x matrices
/// and the contemporaneous matrices A_0, unit lower triangular.
struct GridDraw
{
    std::vector<arma::mat> b;
    std::vector<arma::mat> a0;
};

GridDraw grid_draw(const VarSpec &spec, const arma::mat &a, const arma::uword draw,
                   const arma::uword n_x)
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword n_levels = spec.quantiles.size();
    const arma::uword n_structural = static_cast<arma::uword>(spec.n_structural());
    const arma::uword nparams = a.n_rows / n_levels;
    GridDraw out;
    for (arma::uword j = 0; j < n_levels; j++)
    {
        const arma::vec level = a.col(draw).subvec(j * nparams, (j + 1) * nparams - 1);
        out.b.push_back(n_x > 0 ? arma::mat(arma::reshape(level.head(k * n_x), k, n_x))
                                : arma::mat(k, 0));
        arma::mat a0 = arma::eye<arma::mat>(k, k);
        if (spec.structural && n_structural > 0)
        {
            core::fill_strict_lower_triangle_by_column(a0, level.tail(n_structural));
        }
        out.a0.push_back(a0);
    }
    return out;
}

} // namespace

using core::AldShape;
using core::ald_log_density;
using core::ald_shape;
using core::BvsBlock;
using core::BvsScope;
using core::bvs_sweep;
using core::draw_ald_scale;
using core::draw_ald_weights;
using core::draw_normal_precision;
using core::iid_block;
using core::IidBlock;
using core::stacked_response;
using core::report_flat_selection_prior;

VarNormalAldDraws VarNormalAldSampler::draw_coefficients(const VarNormalAldInput &input,
                                                        Reporter &reporter) const
{
    input.validate();

    // A grid: every level estimated as the single-quantile model it is, one
    // after the other, and the draws stacked level by level. Each level's run
    // is the chain a file with that one quantile would have drawn.
    if (input.spec.uses_quantile_grid())
    {
        const std::vector<double> &levels = input.spec.quantiles;
        const long long per_level = input.spec.draws();
        const long long total = per_level * static_cast<long long>(levels.size());
        VarNormalAldDraws out;
        for (std::size_t j = 0; j < levels.size(); j++)
        {
            VarNormalAldInput one = input;
            one.spec.quantiles.clear();
            one.spec.quantile = levels[j];
            one.spec.h = 0;
            one.spec.forecast_quantile = 0.0;
            one.forecast = ForecastData();
            GridReporter level_reporter(reporter, per_level * static_cast<long long>(j), total);
            const VarNormalAldDraws level = draw_coefficients(one, level_reporter);
            out.a = arma::join_cols(out.a, level.a);
            out.a_lambda = arma::join_cols(out.a_lambda, level.a_lambda);
            out.u_scale = arma::join_cols(out.u_scale, level.u_scale);
        }
        reporter.finish();
        return out;
    }

    const int k = input.spec.k;
    const int iterations = input.spec.iterations;
    const int draws = input.spec.draws();

    const arma::vec y = stacked_response(input.train);
    arma::mat z = input.train.z;

    const int nparams = static_cast<int>(z.n_cols);
    const IidBlock iid = iid_block(input.spec, static_cast<arma::uword>(nparams));
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
    arma::vec a, a_prior_rhs;
    arma::mat a_prior_vinv, a_post_v;

    // Variable selection
    std::optional<BvsBlock> a_bvs;
    arma::mat z_bvs;

    if (use_a)
    {
        // The i.i.d. block, if the model has one: the restricted equations'
        // columns leave the regressors, their rows and columns leave the prior,
        // and the zeros go back when a draw is stored. Every one of these is the
        // identity when nothing is restricted.
        //
        // The prior's precision-weighted mean is formed once here rather than
        // once per draw. It was a constant inside the loop before, so the draws
        // are unchanged; reduced, it is the free part of V mu, which is what the
        // prior conditional on the restricted coefficients being zero asks for.
        a_prior_rhs = iid.elements(input.a_prior.v_inv * input.a_prior.mu);
        a_prior_vinv = iid.block(input.a_prior.v_inv);
        a = iid.elements(input.initial.a);
        z = iid.columns(z);
        out.a = arma::mat(nparams, iterations);

        if (use_varsel)
        {
            out.a_lambda = arma::mat(nparams, iterations);

            if (use_bvs)
            {
                z_bvs = z;
                a_bvs.emplace(input.initial.a_lambda, input.a_varsel_prior);
                report_flat_selection_prior(reporter, input.spec.varsel, "a", input.a_varsel_prior,
                                            input.a_prior.v_inv);
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
            a = draw_normal_precision(a_post_v, a_prior_rhs +
                                                    arma::trans(z) * u_sigma_inv_diag * y_adjusted);

            if (use_bvs)
            {
                // The skew has to be inside the residual here as well. Scoring
                // an inclusion against y - z*theta rather than
                // y - theta*w - z*theta would select for the median whatever
                // quantile was asked for, and at q = 0.5 nothing would show.
                //
                // Against the unmasked regressors, as every other sampler scores
                // them: the sweep hands over candidates already masked. Masking
                // them a second time, with the indicators the sweep is still
                // updating, zeroed the "on" candidate of any coefficient that
                // was out, so its way back in was decided by the prior alone.
                bvs_sweep(*a_bvs, a, BvsScope::element, [&](const arma::vec &theta) {
                    const arma::vec res = y_adjusted - z_bvs * theta;
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
        if (input.spec.keeps(draw))
        {
            const int draw_pos = input.spec.kept_index(draw);

            if (use_a)
            {
                out.a.col(draw_pos) = iid.scatter(a);
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
    if (!input.spec.uses_quantile_grid())
    {
        throw std::invalid_argument(
            "a quantile regression model does not forecast: the h step ahead quantile is not the "
            "quantile of the iterated one step ahead quantiles, so iterating this model forward "
            "would produce a path that cannot be read as a quantile of anything. A grid of "
            "quantiles, /model/quantiles, describes the whole distribution and does forecast");
    }
    input.validate();

    const int k = input.spec.k;
    const int p = input.spec.p;
    const int h = input.spec.h;
    if (h <= 0)
    {
        throw std::invalid_argument("forecast horizon (h) must be positive");
    }
    if (draws.u_scale.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of the asymmetric Laplace scale are missing");
    }

    arma::mat x = input.forecast.x;
    core::require_forecast_regressors(input.spec, x);
    core::require_forecast_horizons(x, h);

    const arma::uword n_levels = input.spec.quantiles.size();
    const arma::vec tau(input.spec.quantiles);
    const bool use_a = x.n_elem > 0 || input.spec.structural;
    if (use_a && (!draws.has_a() || draws.a.n_rows % n_levels != 0))
    {
        throw std::invalid_argument("posterior draws of a are missing, or do not hold one block per "
                                    "level of /model/quantiles");
    }
    if (x.n_elem > 0 && draws.a.n_rows / n_levels !=
                            x.n_cols * static_cast<arma::uword>(k) +
                                static_cast<arma::uword>(input.spec.n_structural()))
    {
        throw std::invalid_argument(
            "forecast regressors and coefficient draws disagree: x has " +
            std::to_string(x.n_cols) + " columns, and each level of a has " +
            std::to_string(draws.a.n_rows / n_levels) + " coefficients");
    }

    // The scenario: a pinned value replaces the draw of that variable in that
    // period, and everything ordered after it and every later period responds.
    arma::mat pinned(k, h);
    pinned.fill(arma::datum::nan);
    const Constraints &c = input.forecast.constraints;
    for (arma::uword e = 0; e < c.row.n_elem; e++)
    {
        pinned(c.variable(e), c.period(e)) = c.value(c.row(e)) / c.weight(e);
    }

    // Every draw of a level is used, so the level is drawn even where a pin
    // makes it unused: a scenario and the same forecast without it then share
    // their random numbers, and their difference is the response to the pins.
    const bool fixed_level = input.spec.forecast_quantile > 0.0;
    const arma::uword n_draws = draws.u_scale.n_cols;
    arma::mat fcst(static_cast<arma::uword>(h * k), n_draws, arma::fill::zeros);
    arma::vec q(n_levels);

    for (arma::uword draw = 0; draw < n_draws; draw++)
    {
        reporter.check_interrupt();
        reporter.progress(static_cast<long long>(draw) + 1, static_cast<long long>(n_draws));

        const GridDraw g = use_a ? grid_draw(input.spec, draws.a, draw, x.n_cols) : GridDraw();
        for (int i = 0; i < h; i++)
        {
            if (i > 0 && p > 0 && x.n_elem > 0)
            {
                core::update_forecast_lags(x, fcst, draw, i, k, p);
            }
            arma::vec y(static_cast<arma::uword>(k), arma::fill::zeros);
            for (int v = 0; v < k; v++)
            {
                for (arma::uword j = 0; j < n_levels; j++)
                {
                    double value = 0.0;
                    if (use_a)
                    {
                        if (x.n_elem > 0)
                        {
                            value = arma::dot(g.b[j].row(v), x.row(i));
                        }
                        // A_0 y = B x + u: each earlier variable enters with minus
                        // its contemporaneous coefficient.
                        for (int w = 0; w < v; w++)
                        {
                            value -= g.a0[j](v, w) * y(w);
                        }
                    }
                    q(j) = value;
                }
                const double level = fixed_level ? input.spec.forecast_quantile : open_uniform();
                const double pin = pinned(v, i);
                const double lower_scale = draws.u_scale(static_cast<arma::uword>(v), draw);
                const double upper_scale =
                    draws.u_scale((n_levels - 1) * static_cast<arma::uword>(k) + v, draw);
                y(v) = std::isnan(pin)
                           ? core::quantile_grid_value(
                                 core::quantile_grid(tau, q, lower_scale, upper_scale), level)
                           : pin;
            }
            fcst.submat(i * k, draw, (i + 1) * k - 1, draw) = y;
        }
    }

    reporter.finish();
    ForecastDraws out;
    out.values = fcst;
    return out;
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

    // A grid: the density its quantiles describe, at every observation. The
    // contemporaneous terms are regressors in z, so one product per level gives
    // every conditional quantile of the sample.
    if (input.spec.uses_quantile_grid())
    {
        const arma::uword n_levels = input.spec.quantiles.size();
        const arma::vec tau(input.spec.quantiles);
        if (use_a && coefficients.a.n_rows != z.n_cols * n_levels)
        {
            throw std::invalid_argument("posterior draws of a do not hold one block per level of "
                                        "/model/quantiles");
        }
        const arma::uword nparams = z.n_cols;
        arma::mat loglik(draws, static_cast<arma::uword>(tt), arma::fill::zeros);
        arma::mat levels(y.n_elem, n_levels, arma::fill::zeros);
        for (arma::uword draw = 0; draw < draws; draw++)
        {
            if (use_a)
            {
                for (arma::uword j = 0; j < n_levels; j++)
                {
                    levels.col(j) =
                        z * coefficients.a.col(draw).subvec(j * nparams, (j + 1) * nparams - 1);
                }
            }
            for (int t = 0; t < tt; t++)
            {
                double total = 0.0;
                for (int v = 0; v < k; v++)
                {
                    const arma::uword r = static_cast<arma::uword>(t * k + v);
                    const double lower_scale = coefficients.u_scale(static_cast<arma::uword>(v), draw);
                    const double upper_scale = coefficients.u_scale(
                        (n_levels - 1) * static_cast<arma::uword>(k) + v, draw);
                    total += core::quantile_grid_log_density(
                        core::quantile_grid(tau, arma::vec(levels.row(r).t()), lower_scale, upper_scale),
                        y(r));
                }
                loglik(draw, static_cast<arma::uword>(t)) = total;
            }
        }
        return loglik;
    }

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
