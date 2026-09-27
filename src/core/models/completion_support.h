// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_COMPLETION_SUPPORT_H
#define BAYESTS_CORE_MODELS_COMPLETION_SUPPORT_H

#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/spec.h"
#include "core/algorithms/constrained_var_path.h"
#include "core/models/constraint_support.h"
#include "core/models/model_support.h"
#include "core/models/predictive_score.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace bayests::core
{

// What a VAR sampler needs to run the data-completion step in its sweep. The
// step itself is draw_constrained_var_path(); these are the pieces on either
// side of it: the regressors a completed panel implies, and the prior over the
// path a draw of the parameters implies.
//
// The lag block is the first k p columns of a regressor row, lag j at columns
// (j - 1) k ... j k - 1, the layout update_forecast_lags() writes. The rest of
// the row -- exogenous and deterministic terms -- is data, and a completed panel
// changes none of it.

/// The compact regressors, tt x n_x, one row per period, read back off the SUR
/// matrix the file holds: z is kron(x, I_k), so x(t, j) is its element at row
/// t k and column j k.
inline arma::mat compact_regressors(const arma::mat &z, const arma::uword k)
{
    const arma::uword tt = z.n_rows / k;
    const arma::uword n_x = z.n_cols / k;
    arma::mat x(tt, n_x);
    for (arma::uword t = 0; t < tt; t++)
    {
        for (arma::uword j = 0; j < n_x; j++)
        {
            x(t, j) = z(t * k, j * k);
        }
    }
    return x;
}

/// The p periods before the sample, k x p and oldest first, as the first row of
/// regressors holds them: y_{-m} is lag m of period zero.
///
/// These are conditioned on, as a VAR conditions on its first p observations,
/// so the host fills them for a series whose presample it did not observe --
/// interpolated, say. Drawing them too would need a prior on the state before
/// the sample, which is a later refinement.
inline arma::mat lag_presample(const arma::mat &x, const arma::uword k, const arma::uword p)
{
    arma::mat presample(k, p);
    for (arma::uword m = 1; m <= p; m++)
    {
        presample.col(p - m) = arma::trans(x.submat(0, (m - 1) * k, 0, m * k - 1));
    }
    return presample;
}

/// Overwrites the lag block of every row whose lag reaches into the sample with
/// the completed panel, k x tt. Rows whose lag reaches before the sample keep
/// what they hold: that is the presample, which is conditioned on.
inline void fill_lags(arma::mat &x, const arma::mat &path, const arma::uword p)
{
    const arma::uword k = path.n_rows;
    for (arma::uword t = 1; t < x.n_rows; t++)
    {
        for (arma::uword j = 1; j <= p && j <= t; j++)
        {
            x.submat(t, (j - 1) * k, t, j * k - 1) = arma::trans(path.col(t - j));
        }
    }
}

/// The Gaussian a draw of a constant-coefficient VAR puts on its own path:
/// `a_matrix` is the draw's k x n_x coefficients, `x` the compact regressors
/// (only the columns past the lag block are read), `covariance` the error
/// covariance.
inline VarPathPrior var_path_prior(const arma::mat &a_matrix, const arma::mat &x,
                                   const arma::mat &presample, const arma::mat &covariance,
                                   const arma::uword p)
{
    // Off the coefficients, not the covariance: a covariance that moves comes as
    // one block per period, tt k rows of them.
    const arma::uword k = a_matrix.n_rows;
    const arma::uword lag = k * p;
    const arma::uword n_x = a_matrix.n_cols;

    VarPathPrior prior;
    prior.presample = presample;
    prior.coefficients = lag > 0 ? arma::mat(a_matrix.cols(0, lag - 1)) : arma::mat(k, 0);
    prior.offset = n_x > lag ? arma::mat(a_matrix.cols(lag, n_x - 1) * x.cols(lag, n_x - 1).t())
                             : arma::mat(k, x.n_rows, arma::fill::zeros);
    prior.covariance = covariance;
    return prior;
}

/// The last p periods of a panel, k x p and oldest first: where a horizon
/// starts. `panel` is k x tt.
inline arma::mat horizon_presample(const arma::mat &panel, const arma::uword p)
{
    if (p == 0)
    {
        return arma::mat(panel.n_rows, 0);
    }
    return panel.cols(panel.n_cols - p, panel.n_cols - 1);
}

/// A constraint set pinning every entry of `realised`, one row per period and
/// one column per variable: a horizon realised whole, in the form a horizon
/// realised in part takes. Its log density under a path's prior is the
/// pointwise log likelihood of those periods.
inline Constraints pin_every_entry(const arma::mat &realised)
{
    const arma::uword periods = realised.n_rows;
    const arma::uword k = realised.n_cols;
    const arma::uword n = periods * k;
    Constraints c;
    c.value.set_size(n);
    c.group.zeros(n);
    c.row = arma::regspace<arma::uvec>(0, n - 1);
    c.period.set_size(n);
    c.variable.set_size(n);
    c.weight.ones(n);
    for (arma::uword t = 0; t < periods; t++)
    {
        for (arma::uword i = 0; i < k; i++)
        {
            const arma::uword r = t * k + i;
            c.value(r) = realised(t, i);
            c.period(r) = t;
            c.variable(r) = i;
        }
    }
    return c;
}

/// How many soft groups a constraint set has: its largest group, soft groups
/// being numbered 1..G with none left out.
inline arma::uword soft_groups(const Constraints &c)
{
    return c.group.n_elem > 0 ? c.group.max() : 0;
}

/// The error variances the completion step is handed, one per soft group, from
/// their precisions.
inline arma::vec soft_variances(const arma::vec &precision)
{
    return precision.n_elem > 0 ? arma::vec(1.0 / precision) : arma::vec();
}

/// One Gibbs draw of each soft group's error precision given the completed
/// panel, k x tt: a gamma prior is conjugate, so group g's posterior has shape
/// shape_g + n_g / 2 and rate rate_g + (sum of its rows' squared errors) / 2,
/// the error of a row being its value less what the panel says of it. Hard rows
/// hold exactly and contribute nothing.
///
/// Consumes one gamma draw per group and none where no row is soft, so a set of
/// hard rows moves no random number.
inline void draw_soft_precisions(arma::vec &precision, const Constraints &c, const arma::mat &path,
                                 const GammaPrior &prior)
{
    const arma::uword groups = soft_groups(c);
    if (groups == 0)
    {
        return;
    }

    arma::vec fitted(c.value.n_elem, arma::fill::zeros);
    for (arma::uword e = 0; e < c.row.n_elem; e++)
    {
        fitted(c.row(e)) += c.weight(e) * path(c.variable(e), c.period(e));
    }
    arma::vec sse(groups, arma::fill::zeros);
    arma::vec count(groups, arma::fill::zeros);
    for (arma::uword r = 0; r < c.value.n_elem; r++)
    {
        if (c.group(r) > 0)
        {
            const double error = c.value(r) - fitted(r);
            sse(c.group(r) - 1) += error * error;
            count(c.group(r) - 1) += 1.0;
        }
    }
    for (arma::uword g = 0; g < groups; g++)
    {
        precision(g) = arma::randg<double>(arma::distr_param(
            prior.shape(g) + 0.5 * count(g), 1.0 / (prior.rate(g) + 0.5 * sse(g))));
    }
}

/// Refuses soft rows in the training set without a prior and a start for their
/// precisions: one shape, one rate and one initial precision per group, each
/// positive and finite. Nothing is asked of a set without soft rows.
inline void require_soft_prior(const Constraints &c, const GammaPrior &prior,
                               const arma::vec &initial, const char *algorithm)
{
    const arma::uword groups = soft_groups(c);
    if (groups == 0)
    {
        return;
    }
    const auto positive = [&](const arma::vec &v, const std::string &what) {
        if (v.n_elem != groups)
        {
            throw std::invalid_argument(
                std::string(algorithm) + ": /data/train/constraints has " + std::to_string(groups) +
                " soft groups, so " + what + " needs one value per group, got " +
                std::to_string(v.n_elem));
        }
        require_finite(v, what);
        if (!(v.min() > 0.0))
        {
            throw std::invalid_argument(what + " must be positive in every group");
        }
    };
    positive(prior.shape, "/priors/constraints/shape");
    positive(prior.rate, "/priors/constraints/rate");
    positive(initial, "/initial/constraints_inv");
}

/// var_path_prior() from a draw of the coefficients as the samplers carry it:
/// one column vec(A) of the k x n_x coefficient matrix for a constant VAR, or
/// one such column per period, nparams x tt, for a time-varying one. The
/// covariance is k x k, or one block per period stacked.
inline VarPathPrior var_path_prior_from(const arma::mat &a, const arma::uword k, const arma::mat &x,
                                        const arma::mat &presample, const arma::mat &covariance,
                                        const arma::uword p)
{
    const arma::uword n_x = x.n_cols;
    if (a.n_cols <= 1)
    {
        return var_path_prior(n_x > 0 ? arma::mat(arma::reshape(a, k, n_x)) : arma::mat(k, 0), x,
                              presample, covariance, p);
    }

    const arma::uword tt = x.n_rows;
    const arma::uword lag = k * p;
    VarPathPrior prior;
    prior.presample = presample;
    prior.coefficients = arma::mat(tt * k, lag);
    prior.offset = arma::mat(k, tt, arma::fill::zeros);
    prior.covariance = covariance;
    for (arma::uword t = 0; t < tt; t++)
    {
        const arma::mat a_t = arma::reshape(a.col(t), k, n_x);
        if (lag > 0)
        {
            prior.coefficients.rows(t * k, t * k + k - 1) = a_t.cols(0, lag - 1);
        }
        if (n_x > lag)
        {
            prior.offset.col(t) = a_t.cols(lag, n_x - 1) * arma::trans(x.submat(t, lag, t, n_x - 1));
        }
    }
    return prior;
}

/// The covariance a precision stands for, the completion step taking the one
/// and the samplers carrying the other. `algorithm` names the sampler in the
/// message.
inline arma::mat covariance_of(const arma::mat &precision, const char *algorithm)
{
    arma::mat covariance;
    if (!arma::inv_sympd(covariance, precision))
    {
        throw std::runtime_error(std::string(algorithm) +
                                 ": the error precision is not positive definite, so the panel "
                                 "cannot be completed");
    }
    return covariance;
}

/// The data-completion step in a VAR sampler's sweep, and what it keeps from
/// one sweep to the next: the compact regressors whose lag block it rewrites,
/// the presample it conditions on, and the error precision of each soft group.
///
/// A sampler constructs one before its first sweep and, when active(), calls
/// complete() at the start of every sweep with its current coefficients and
/// error covariance; what it gets back is the panel to run the rest of the
/// sweep on, and regressors() the SUR matrix that goes with it. Inactive --
/// no constraints -- it does nothing and consumes no random number, which is
/// what keeps a complete panel's draws where they were.
class PanelCompletion
{
public:
    PanelCompletion(const VarSpec &spec, const TrainData &train, const bool use_a,
                    const GammaPrior &soft_prior, const arma::vec &initial_soft)
        : constraints_(train.constraints), soft_prior_(soft_prior),
          k_(static_cast<arma::uword>(spec.k)), p_(static_cast<arma::uword>(spec.p)),
          use_a_(use_a), active_(!train.constraints.empty())
    {
        if (!active_)
        {
            return;
        }
        const arma::uword tt = train.y.n_elem / k_;
        x_ = use_a_ ? compact_regressors(train.z, k_) : arma::mat(tt, 0);
        presample_ = lag_presample(x_, k_, p_);
        soft_ = core::soft_groups(train.constraints);
        precision_ = soft_ > 0 ? initial_soft : arma::vec();
    }

    bool active() const { return active_; }
    arma::uword soft_groups() const { return soft_; }

    /// The compact regressors and presample, for a stage after the chain -- a
    /// log likelihood -- that needs the prior over the path again.
    const arma::mat &regressors_compact() const { return x_; }
    const arma::mat &presample() const { return presample_; }

    /// Draws the panel given the sweep's coefficients (see var_path_prior_from())
    /// and error covariance, then each soft group's precision given the panel,
    /// and rebuilds the lag block of the regressors from it. Returns the panel,
    /// k x tt.
    arma::mat complete(const arma::mat &a, const arma::mat &covariance)
    {
        const arma::mat path = draw_constrained_var_path(
            var_path_prior_from(a, k_, x_, presample_, covariance, p_), constraints_,
            soft_variances(precision_));
        draw_soft_precisions(precision_, constraints_, path, soft_prior_);
        if (use_a_)
        {
            fill_lags(x_, path, p_);
        }
        return path;
    }

    /// The SUR regressors of the last completed panel, (tt k) x nparams.
    arma::mat regressors() const
    {
        return x_.n_cols > 0 ? arma::mat(arma::kron(x_, arma::eye<arma::mat>(k_, k_))) : arma::mat();
    }

    /// Sizes the two outputs a chain keeps of this: the completed panel and the
    /// soft precisions, one column per kept draw. Left empty when inactive.
    void allocate(arma::mat &y_draws, arma::mat &soft_draws, const arma::uword iterations,
                  const arma::uword tt) const
    {
        if (active_)
        {
            y_draws = arma::mat(k_ * tt, iterations);
        }
        if (soft_ > 0)
        {
            soft_draws = arma::mat(soft_, iterations);
        }
    }

    void store(arma::mat &y_draws, arma::mat &soft_draws, const arma::uword position,
               const arma::vec &y) const
    {
        if (active_)
        {
            y_draws.col(position) = y;
        }
        if (soft_ > 0)
        {
            soft_draws.col(position) = precision_;
        }
    }

private:
    const Constraints &constraints_;
    const GammaPrior &soft_prior_;
    arma::uword k_;
    arma::uword p_;
    bool use_a_;
    bool active_;
    arma::mat x_;
    arma::mat presample_;
    arma::uword soft_ = 0;
    arma::vec precision_;
};

/// The pointwise log likelihood of a panel not observed whole, one draw: the
/// density of what the constraints observed with the rest integrated out,
/// period by period, given the draw's coefficients `a` (see
/// var_path_prior_from()), error covariance and soft precisions -- empty where
/// no row is soft.
inline arma::rowvec observed_log_likelihood(const Constraints &constraints, const arma::mat &x,
                                            const arma::mat &presample, const arma::mat &a,
                                            const arma::mat &covariance,
                                            const arma::vec &soft_precision, const arma::uword k,
                                            const arma::uword p)
{
    return arma::trans(constrained_var_path_log_density(
        var_path_prior_from(a, k, x, presample, covariance, p), constraints,
        soft_variances(soft_precision)));
}

/// Refuses a chain's output that a panel with soft rows cannot be read without:
/// the precision of each soft group, one row per group and one column per draw.
inline void require_soft_draws(const Constraints &constraints, const arma::mat &soft_draws,
                               const arma::uword draws)
{
    const arma::uword soft = soft_groups(constraints);
    if (soft > 0 && (soft_draws.n_rows != soft || soft_draws.n_cols != draws))
    {
        throw std::invalid_argument(
            "/data/train/constraints has " + std::to_string(soft) + " soft groups, and the "
            "density of what they observed needs their precisions: "
            "/posterior/constraints_inv/coeffs is missing or not one row per group and one "
            "column per draw");
    }
}

/// The lags of a forecast's first horizons that reach back into the sample,
/// from a completed panel, k x tt: horizon i's lag j for every j > i. The lags
/// inside the horizon are the forecast's own, written as it unfolds.
inline void start_from_panel(arma::mat &x, const arma::mat &panel, const arma::uword p)
{
    const arma::uword k = panel.n_rows;
    const arma::uword tt = panel.n_cols;
    for (arma::uword i = 0; i < std::min<arma::uword>(x.n_rows, p); i++)
    {
        for (arma::uword j = i + 1; j <= p; j++)
        {
            x.submat(i, (j - 1) * k, i, j * k - 1) = arma::trans(panel.col(tt + i - j));
        }
    }
}

// What the constant-coefficient VARs share on the other side of the chain --
// the log likelihood, the forecast's start, a scenario and a score -- written
// once against their draws: `a` one column vec(A) per draw, `u_sigma_inv` one
// column of the whole k x k precision per draw, `y` the completed panels and
// `constraints_inv` the soft precisions. VarNormalWishart and VarNormalGamma
// both carry exactly that.

/// The panel a forecast from draw `draw` starts from, k x tt: the draw's own
/// completion where the sample was not observed whole, the sample otherwise.
template <typename Input, typename Draws>
arma::mat starting_panel(const Input &input, const Draws &draws, const arma::uword draw,
                         const arma::uword k)
{
    const arma::vec y = stacked_response(input.train);
    const arma::uword tt = y.n_elem / k;
    return input.train.constraints.empty() ? arma::mat(arma::reshape(y, k, tt))
                                           : arma::mat(arma::reshape(draws.y.col(draw), k, tt));
}

/// The error covariance of every period from a stack of k x k precision
/// blocks, one per period: the stack the time-varying samplers carry.
inline arma::mat stacked_covariance(const arma::mat &blocks, const arma::uword k,
                                    const char *algorithm)
{
    const arma::uword tt = blocks.n_rows / k;
    arma::mat out(tt * k, k);
    for (arma::uword t = 0; t < tt; t++)
    {
        out.rows(t * k, t * k + k - 1) =
            covariance_of(arma::mat(blocks.rows(t * k, t * k + k - 1)), algorithm);
    }
    return out;
}

/// The same from a block-diagonal precision, (tt k) x (tt k), the form the
/// stochastic volatility sampler carries it in.
inline arma::mat stacked_covariance(const arma::sp_mat &block_diagonal, const arma::uword k,
                                    const char *algorithm)
{
    const arma::uword tt = block_diagonal.n_rows / k;
    arma::mat blocks(tt * k, k);
    for (arma::uword t = 0; t < tt; t++)
    {
        blocks.rows(t * k, t * k + k - 1) =
            arma::mat(block_diagonal.submat(t * k, t * k, t * k + k - 1, t * k + k - 1));
    }
    return stacked_covariance(blocks, k, algorithm);
}

/// Draw `draw`'s error covariance, from the precision the draws carry: one
/// k x k block, or one per period, which is what precision_stride() tells
/// apart.
template <typename Draws>
arma::mat draw_covariance(const Draws &draws, const arma::uword draw, const arma::uword k,
                          const char *algorithm, const arma::uword tt = 1)
{
    const arma::uword stride =
        precision_stride(draws.u_sigma_inv, static_cast<int>(k), static_cast<int>(tt));
    if (stride == 0)
    {
        return covariance_of(arma::reshape(draws.u_sigma_inv.col(draw), k, k), algorithm);
    }
    arma::mat blocks(tt * k, k);
    for (arma::uword t = 0; t < tt; t++)
    {
        blocks.rows(t * k, t * k + k - 1) =
            arma::reshape(draws.u_sigma_inv.submat(t * stride, draw, t * stride + k * k - 1, draw),
                          k, k);
    }
    return stacked_covariance(blocks, k, algorithm);
}

/// Draw `draw`'s coefficients as var_path_prior_from() takes them, or none: one
/// column for a constant VAR, and a path, nparams x tt, where the column holds
/// one set per period.
template <typename Draws>
arma::mat draw_coefficients_column(const Draws &draws, const arma::uword draw, const bool use_a,
                                   const arma::uword nparams = 0, const arma::uword tt = 1)
{
    if (!use_a)
    {
        return arma::mat();
    }
    if (tt > 1 && nparams > 0 && draws.a.n_rows == nparams * tt)
    {
        return arma::reshape(draws.a.col(draw), nparams, tt);
    }
    return arma::mat(draws.a.col(draw));
}

/// Refuses a score where the sample was not observed whole, or the horizon was
/// realised in part, for a model that does not score from a completed panel
/// yet. A score runs without validate() in between, so this is where it is
/// refused rather than computed from the placeholders.
template <typename Input>
void require_no_panel_score(const Input &input, const char *algorithm)
{
    if (!input.train.constraints.empty() || !input.test.constraints.empty())
    {
        throw std::invalid_argument(
            std::string(algorithm) +
            " does not score a forecast from a panel not observed whole, or against a horizon "
            "realised in part, yet");
    }
}

/// The pointwise log likelihood of a panel not observed whole, every draw: the
/// density of what was observed, with the rest integrated out.
template <typename Input, typename Draws>
arma::mat observed_log_likelihood(const Input &input, const Draws &draws, const char *algorithm)
{
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword p = static_cast<arma::uword>(input.spec.p);
    const bool use_a = input.train.z.n_cols > 0;
    const arma::uword tt = input.train.y.n_elem / k;
    const arma::mat x = use_a ? compact_regressors(input.train.z, k) : arma::mat(tt, 0);
    const arma::mat presample = lag_presample(x, k, p);
    const arma::uword n_draws = draws.iterations();
    require_soft_draws(input.train.constraints, draws.constraints_inv, n_draws);
    const bool soft = soft_groups(input.train.constraints) > 0;

    const arma::uword nparams = input.train.z.n_cols;
    arma::mat loglik(n_draws, tt);
    for (arma::uword draw = 0; draw < n_draws; draw++)
    {
        loglik.row(draw) = observed_log_likelihood(
            input.train.constraints, x, presample,
            draw_coefficients_column(draws, draw, use_a, nparams, tt),
            draw_covariance(draws, draw, k, algorithm, tt),
            soft ? arma::vec(draws.constraints_inv.col(draw)) : arma::vec(), k, p);
    }
    return loglik;
}

/// Refuses a forecast from a panel not observed whole without the completed
/// panels it starts from.
template <typename Input, typename Draws>
void require_completed_panels(const Input &input, const Draws &draws, const bool use_a)
{
    if (input.train.constraints.empty() || !use_a || input.spec.p <= 0)
    {
        return;
    }
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword tt = input.train.y.n_elem / k;
    if (draws.y.n_rows != tt * k || draws.y.n_cols != draws.iterations())
    {
        throw std::invalid_argument(
            "the panel was not observed whole, so a forecast starts from each draw's completed "
            "panel, and /posterior/y/coeffs is missing or not one column of k tt per draw");
    }
}

/// One draw of a forecast conditioned on /data/forecast/constraints: the path
/// the draw's VAR puts on the horizon, starting where its panel ends, drawn by
/// the completion step. k x h.
template <typename Input, typename Draws>
arma::mat conditioned_forecast(const Input &input, const Draws &draws, const arma::uword draw,
                               const bool use_a, const char *algorithm)
{
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword p = static_cast<arma::uword>(input.spec.p);
    const arma::mat x = use_a ? input.forecast.x : arma::mat(static_cast<arma::uword>(input.spec.h), 0);
    return draw_constrained_var_path(
        var_path_prior_from(draw_coefficients_column(draws, draw, use_a), k, x,
                            horizon_presample(starting_panel(input, draws, draw, k), p),
                            draw_covariance(draws, draw, k, algorithm), p),
        input.forecast.constraints, arma::vec());
}

/// The score where the sample was not observed whole or the horizon was
/// realised in part: the completion step's log density of what the horizon
/// realised, from where each draw's panel ends, draws x periods. `periods` is
/// what scored_horizons() found; the realised rows are checked here, a score
/// running from draws without validate() in between.
template <typename Input, typename Draws>
arma::mat score_from_panel(const Input &input, const Draws &draws, const arma::uword periods,
                           const char *algorithm)
{
    const arma::uword k = static_cast<arma::uword>(input.spec.k);
    const arma::uword p = static_cast<arma::uword>(input.spec.p);
    const arma::mat realised = input.test.y.head_rows(periods);

    Constraints observed = pin_every_entry(realised);
    if (!input.test.constraints.empty())
    {
        validate_constraints(input.test.constraints, arma::trans(realised),
                             "/data/test/constraints");
        observed = input.test.constraints;
    }

    const bool use_a = input.train.z.n_cols > 0;
    if (use_a && draws.a.n_elem == 0)
    {
        throw std::invalid_argument("the model has regressors but posterior draws of a are missing");
    }
    if (!input.train.constraints.empty() && p > 0 && draws.y.n_elem == 0)
    {
        throw std::invalid_argument(
            "the sample was not observed whole, so a score starts from each draw's completed "
            "panel, and /posterior/y/coeffs is missing");
    }
    const arma::mat x = use_a ? arma::mat(input.forecast.x.head_rows(periods)) : arma::mat(periods, 0);

    const arma::uword n_draws = draws.iterations();
    arma::mat score(n_draws, periods);
    for (arma::uword draw = 0; draw < n_draws; draw++)
    {
        score.row(draw) = arma::trans(constrained_var_path_log_density(
            var_path_prior_from(draw_coefficients_column(draws, draw, use_a), k, x,
                                horizon_presample(starting_panel(input, draws, draw, k), p),
                                draw_covariance(draws, draw, k, algorithm), p),
            observed, arma::vec()));
    }
    return score;
}

/// The combinations a VAR sampler reading constraints refuses, beside the
/// checks require_supported_constraints() has already made on the sets
/// themselves. `forecast` carries the scenario, `test` what the horizon
/// realised.
///
/// Each is refused rather than approximated. Soft rows on the horizon would
/// need a variance nothing draws or gives. A score of a conditional
/// forecast is refused for what it would mean: a score is the density of what
/// the horizon realised under the model's own forecast, and a scenario replaces
/// that forecast with another. A structural form or an i.i.d. block rearranges
/// the coefficients the lag block is read from.
///
/// `horizon` says whether the sampler scores from a completed panel and
/// conditions a forecast on a scenario. Where it does not -- a model whose
/// forecast carries states forward, until that is done for it -- a score beside
/// training constraints, constraints on what the horizon realised and a
/// scenario are all refused by name.
inline void require_completion_spec(const VarSpec &spec, const TrainData &train,
                                    const ForecastData &forecast, const TestData &test,
                                    const char *algorithm, const bool horizon = true)
{
    if (train.constraints.empty() && test.constraints.empty() && forecast.constraints.empty())
    {
        return;
    }

    const std::string name(algorithm);
    if (!horizon && (!test.constraints.empty() || !forecast.constraints.empty() ||
                     (!train.constraints.empty() && test.y.n_elem > 0)))
    {
        throw std::invalid_argument(
            name + " does not score or condition a forecast from a panel not observed whole "
                   "yet: drop /data/test and /data/forecast/constraints, or estimate without "
                   "/data/train/constraints");
    }

    // Soft rows in the sample have their variances drawn. On the horizon nothing
    // would inform one: a realised value or a scenario with an error would need
    // the error's variance given, which the file has no place for yet.
    const auto require_hard = [&](const Constraints &c, const char *where) {
        if (soft_groups(c) > 0)
        {
            throw std::invalid_argument(
                name + " reads hard constraints only in " + where + ": its group is above zero "
                "somewhere, and nothing draws or gives the variance of a soft row there");
        }
    };
    require_hard(test.constraints, "/data/test/constraints");
    require_hard(forecast.constraints, "/data/forecast/constraints");

    if (!forecast.constraints.empty() && (test.y.n_elem > 0 || !test.constraints.empty()))
    {
        throw std::invalid_argument(
            name + " does not score a conditional forecast: a score is the density of what the "
                   "horizon realised under the model's own forecast, and "
                   "/data/forecast/constraints replaces that forecast with a scenario. Score "
                   "without the scenario, or condition without /data/test");
    }
    if (spec.structural || spec.n_iid > 0)
    {
        throw std::invalid_argument(
            name + " does not combine constraints with a structural form or n_iid: both "
                   "rearrange the coefficients the lags are read from");
    }

    const arma::uword k = spec.k > 0 ? static_cast<arma::uword>(spec.k) : 0;
    const arma::uword want = k * static_cast<arma::uword>(spec.n_x());
    if (train.z.n_cols != want && !(train.z.n_cols == 0 && spec.p == 0))
    {
        throw std::invalid_argument(
            name + " rebuilds the lags of a panel it completes from /data/train/z, so z must "
                   "have k n_x = " + std::to_string(want) + " columns, the lag block first; got " +
            std::to_string(train.z.n_cols));
    }
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_COMPLETION_SUPPORT_H
