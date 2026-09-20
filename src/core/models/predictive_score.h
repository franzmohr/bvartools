// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_PREDICTIVE_SCORE_H
#define BAYESTS_CORE_MODELS_PREDICTIVE_SCORE_H

#include "bayests/arma.h"
#include "bayests/data.h"
#include "bayests/spec.h"

#include "core/algorithms/triangular_packing.h"
#include "core/models/forecast_states.h"

#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Scoring a forecast against what the horizon realised.
///
/// The score of horizon i is the log density of the realised observation under
/// each draw, conditional on the realised observations before it -- not on the
/// path the forecast simulated. That is the choice the whole shape rests on:
/// the log of the mean of exp() over draws is then the one step ahead
/// predictive density given everything known up to that period, and the sum of
/// those over the horizons is log p(y*_{T+1..T+h} | data), the joint. Scoring
/// against a simulated history instead would give the marginal density of each
/// horizon on its own, which is a defensible quantity and a different one --
/// marginals do not sum to a joint, so a column of them could not be added up.
///
/// The consequence that makes this cheap: with the history realised rather than
/// simulated, the regressors of the scored periods do not depend on the draw.
/// They are built once, and the score is then the model's own pointwise log
/// likelihood over those periods -- the same expression, the same code, a
/// different sample. Nothing new is written down per model, which is what keeps
/// one density per algorithm rather than two that can drift apart.

/// How many horizons a file can be scored over: the periods the realised values
/// cover, which read_test_observations() has already refused to let exceed `h`.
///
/// Throws where there is nothing to score. Callers check for realised values
/// before asking, so reaching here without them is a caller's mistake rather
/// than a file's.
inline arma::uword scored_horizons(const arma::mat &realised, const VarSpec &spec)
{
    if (realised.n_rows == 0)
    {
        throw std::invalid_argument(
            "scoring a forecast needs the observations the horizon realised, and /data/test/y "
            "holds none");
    }
    if (spec.h <= 0)
    {
        throw std::invalid_argument(
            "scoring a forecast needs a horizon, and this model asks for none");
    }
    return realised.n_rows;
}

/// Refuses a model this cannot be the density of, saying which part is missing
/// rather than producing a number of the right size and the wrong meaning.
///
/// A structural model is refused, for a reason that does not go away with more
/// code here: its regressors include the contemporaneous observations, so the
/// realised row is not built from lags alone, and its density carries the
/// Jacobian of A_0 besides.
inline void require_scorable(const VarSpec &spec, const std::string &algorithm)
{
    if (spec.structural)
    {
        throw std::invalid_argument(
            algorithm +
            ": a structural model cannot be scored. Its regressors include the contemporaneous "
            "observations and its density carries the Jacobian of A_0, neither of which the "
            "pointwise log likelihood over the scored periods accounts for");
    }
}

/// Each draw's random walk carried forward, one step per scored period, in the
/// layout a pointwise log likelihood reads a path in: `n_state * periods` rows
/// with the periods stacked within a column, one column per draw.
///
/// `last` is the state at the end of the estimation sample, `n_state` per draw,
/// which is what the forecast readers hand over. The step comes before the
/// period it belongs to, exactly as it does in a forecast: period T+1 is
/// generated under the state the first step arrives at, not under the sample's
/// last.
///
/// Under ForecastStates::hold nothing steps and the sample's last state is
/// repeated. That scores the model whose drift stops where the sample does,
/// which is a different model from the one estimated and is worth asking for
/// only deliberately -- see ForecastStates.
///
/// One path per draw is one sample of the state, which is all the density
/// needs: averaging exp() over draws integrates the state out along with
/// everything else the posterior carries. It does mean the score is drawn
/// rather than computed, so two runs of `simulate` give two answers, as two
/// runs of a forecast do; `/model/seed` is what makes either repeatable.
inline arma::mat carry_state_forward(const arma::mat &last, const arma::mat &sigma,
                                     const arma::mat &mask, const arma::uword periods,
                                     const bool simulate, const std::string &name)
{
    const arma::uword n_state = last.n_rows;
    const arma::uword draws = last.n_cols;

    if (simulate)
    {
        require_state_variances(sigma, n_state, draws, name);
        require_state_mask(mask, n_state, draws, name);
    }

    arma::mat path(n_state * periods, draws);
    arma::vec state, step_sigma, step_mask;

    for (arma::uword draw = 0; draw < draws; draw++)
    {
        state = last.col(draw);
        if (simulate)
        {
            step_sigma = sigma.col(draw);
            step_mask = mask.n_elem > 0 ? arma::vec(mask.col(draw)) : arma::vec();
        }
        for (arma::uword i = 0; i < periods; i++)
        {
            if (simulate)
            {
                step_random_walk(state, step_sigma, step_mask);
            }
            path.submat(i * n_state, draw, (i + 1) * n_state - 1, draw) = state;
        }
    }

    return path;
}

/// The same for a stochastic volatility, which is a random walk in the log
/// variances rather than in the precision itself.
///
/// `last_omega_inv` is the diagonal of the precision at the end of the sample,
/// `k` per draw, so the state stepped is `-log()` of it and what comes back is
/// `exp(-h)` again: the precision diagonal, period by period, in the layout
/// `u_omega_inv` is stored in.
inline arma::mat carry_log_volatility_forward(const arma::mat &last_omega_inv,
                                              const arma::mat &h_sigma,
                                              const arma::uword periods, const bool simulate)
{
    const arma::mat log_variance = -arma::log(last_omega_inv);
    const arma::mat path = carry_state_forward(log_variance, h_sigma, arma::mat(), periods,
                                               simulate, "the log-volatilities");
    return arma::exp(-path);
}

/// Psi carried forward. Its free elements are the random walk, so they are
/// packed out of the sample's last Psi, stepped like any other state, and a
/// lower triangular Psi with a unit diagonal is rebuilt at every period.
/// Returns `k * k * periods` rows, one vectorised Psi per period.
///
/// A selection mask arrives as a whole Psi and is packed the same way, which is
/// the order psi_sigma counts its variances in -- see pack_strict_lower_triangle().
inline arma::mat carry_psi_forward(const arma::mat &last_psi, const arma::mat &psi_sigma,
                                   const arma::mat &psi_lambda, const arma::uword periods,
                                   const arma::uword k, const bool simulate)
{
    const arma::uword kk = k * k;
    const arma::uword draws = last_psi.n_cols;
    const arma::uword n_free = k * (k - 1) / 2;

    arma::mat packed(n_free, draws);
    arma::mat mask;
    if (psi_lambda.n_elem > 0)
    {
        mask.set_size(n_free, draws);
    }
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        packed.col(draw) = pack_strict_lower_triangle(arma::reshape(last_psi.col(draw), k, k));
        if (psi_lambda.n_elem > 0)
        {
            mask.col(draw) = pack_strict_lower_triangle(arma::reshape(psi_lambda.col(draw), k, k));
        }
    }

    const arma::mat path = carry_state_forward(packed, psi_sigma, mask, periods, simulate,
                                               "the covariance block");

    arma::mat whole(kk * periods, draws);
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (arma::uword i = 0; i < periods; i++)
        {
            arma::mat psi = arma::eye<arma::mat>(k, k);
            // Materialised rather than passed as a subview: the filler indexes
            // it with subvec(), which a subview does not have.
            const arma::vec free =
                path.submat(i * n_free, draw, (i + 1) * n_free - 1, draw);
            fill_strict_lower_triangle(psi, free);
            whole.submat(i * kk, draw, (i + 1) * kk - 1, draw) = arma::vectorise(psi);
        }
    }

    return whole;
}

/// The precision of each scored period, formed the way the samplers form it:
/// Psi' Omega^-1 Psi, with Omega^-1 the diagonal the volatility gives.
///
/// Both arguments are paths of `periods`, `k` and `k * k` rows per period; a
/// block that stands still is repeated by the caller, which keeps this one
/// expression rather than a stride per argument. An empty `psi_path` is a model
/// with no covariance block, where Psi is the identity and the precision is the
/// diagonal itself.
inline arma::mat precision_path(const arma::mat &omega_inv_path, const arma::mat &psi_path,
                                const arma::uword k, const arma::uword periods)
{
    const arma::uword kk = k * k;
    const arma::uword draws = omega_inv_path.n_cols;
    const bool use_psi = psi_path.n_elem > 0;

    arma::mat precision(kk * periods, draws);
    for (arma::uword draw = 0; draw < draws; draw++)
    {
        for (arma::uword i = 0; i < periods; i++)
        {
            arma::mat period = arma::diagmat(
                omega_inv_path.submat(i * k, draw, (i + 1) * k - 1, draw));
            if (use_psi)
            {
                const arma::mat psi =
                    arma::reshape(psi_path.submat(i * kk, draw, (i + 1) * kk - 1, draw), k, k);
                period = arma::trans(psi) * period * psi;
            }
            precision.submat(i * kk, draw, (i + 1) * kk - 1, draw) = arma::vectorise(period);
        }
    }

    return precision;
}

/// The regressors of the scored periods, with the lagged endogenous blocks
/// taken from what was realised.
///
/// `forecast_x` is `/data/forecast/x`, the compact layout of one period per row
/// that a forecast is driven by, and it comes in with its deterministic and
/// unmodelled columns filled and its lag blocks holding whatever the caller put
/// there. A forecast overwrites those blocks as it simulates; this overwrites
/// them from `realised` instead, and the two do it by the same rule, so a lag
/// reaching back before the first scored period is left exactly as the caller
/// supplied it -- that value is the end of the estimation sample and is already
/// right. See update_forecast_lags(), which this mirrors.
///
/// Rows beyond the scored periods are dropped: a file may forecast further than
/// it realised, and the horizons it did not realise cannot be scored.
inline arma::mat realised_regressors(const arma::mat &forecast_x, const arma::mat &realised,
                                     const int k, const int p)
{
    const arma::uword periods = realised.n_rows;
    if (forecast_x.n_rows < periods)
    {
        throw std::invalid_argument(
            "/data/forecast/x has " + std::to_string(forecast_x.n_rows) +
            " horizons but /data/test/y realised " + std::to_string(periods) +
            "; the regressors of a period that is scored have to be there");
    }
    if (realised.n_cols != static_cast<arma::uword>(k))
    {
        throw std::invalid_argument(
            "the realised values must have k = " + std::to_string(k) + " columns, got " +
            std::to_string(realised.n_cols));
    }

    arma::mat x = forecast_x.head_rows(periods);
    if (p <= 0)
    {
        return x;
    }

    const arma::uword width = static_cast<arma::uword>(k);
    for (arma::uword i = 1; i < periods; i++)
    {
        const arma::uword filled = i < static_cast<arma::uword>(p) ? i : static_cast<arma::uword>(p);
        for (arma::uword j = 1; j <= filled; j++)
        {
            x.submat(i, (j - 1) * width, i, j * width - 1) = realised.row(i - j);
        }
    }

    return x;
}

/// The SUR spelling of those regressors, which is what a sampler's pointwise
/// log likelihood reads: kron(x, I_k), one block row per period.
inline arma::mat sur_regressors(const arma::mat &x, const int k)
{
    if (x.n_cols == 0)
    {
        return arma::mat();
    }
    return arma::kron(x, arma::eye<arma::mat>(k, k));
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_PREDICTIVE_SCORE_H
