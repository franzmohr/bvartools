// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_STEADY_STATE_SUPPORT_H
#define BAYESTS_CORE_MODELS_STEADY_STATE_SUPPORT_H

#include "bayests/priors.h"
#include "bayests/spec.h"
#include "core/models/completion_support.h"
#include "core/models/model_support.h"
#include "core/models/shrinkage_support.h"

#include <functional>
#include <stdexcept>
#include <string>

namespace bayests::core
{

// The steady-state prior of Villani (2009): the VAR in its mean-adjusted form,
//
//     y_t - mu = A_1 (y_{t-1} - mu) + ... + A_p (y_{t-p} - mu) + u_t,
//
// with the prior on mu, the unconditional mean, rather than on the intercept
// (I - sum_j A_j) mu. A prior on where the series settle is one a forecaster can
// state; one on the intercept is not, the intercept's meaning moving with every
// A_j. Two Gibbs blocks take the place of the coefficient draw: the lags given
// mu, a regression of the demeaned series on its demeaned lags, and mu given
// the lags, a regression of y_t - sum_j A_j y_{t-j} on I - sum_j A_j. The
// intercept they imply goes back into `a`, so a forecast, a log likelihood and
// a score read a draw as they read any other.
//
// For a model whose only deterministic term is an intercept and that has no
// exogenous variables: the regressor row is [y_{t-1}' ... y_{t-p}' 1].

/// Which positions of `a` -- vec of the k x (k p + 1) coefficient matrix -- are
/// lags, and which the intercept: the last k.
inline arma::uvec steady_state_lag_positions(const arma::uword k, const arma::uword p)
{
    return arma::regspace<arma::uvec>(0, k * k * p - 1);
}

/// Refuses the steady-state prior on an algorithm that does not read it.
inline void require_supported_steady_state(const VarSpec &spec, const bool supported,
                                           const char *algorithm)
{
    if (spec.steady_state && !supported)
    {
        throw std::invalid_argument(
            std::string(algorithm) +
            " does not read /model/steady_state: only VarNormalWishart, VarNormalGamma and "
            "VarNormalStochvol put their prior on the unconditional mean");
    }
}

/// Refuses the steady-state prior where it is not read or cannot be: a model
/// with regressors other than p lags and an intercept, a prior mean and
/// precision for mu of the wrong shape or not positive definite, a start of the
/// wrong length, a base prior that couples the lags with the intercept, and the
/// combinations that rearrange the coefficient draw.
inline void require_steady_state(const VarSpec &spec, const bool supported, const TrainData &train,
                                 const NormalPrior &a_prior, const NormalPrior &mu_prior,
                                 const arma::vec &initial_mu, const char *algorithm)
{
    if (!spec.steady_state)
    {
        return;
    }
    const std::string name(algorithm);
    if (!supported)
    {
        throw std::invalid_argument(
            name + " does not read /model/steady_state: only VarNormalWishart, VarNormalGamma and "
                   "VarNormalStochvol put their prior on the unconditional mean");
    }
    if (spec.uses_varsel() || spec.structural || spec.n_iid > 0 ||
        spec.shrinkage != Shrinkage::none)
    {
        throw std::invalid_argument(
            name + " does not combine /model/steady_state with variable selection, a structural "
                   "form, n_iid or /model/shrinkage");
    }
    if (spec.p < 1 || spec.n != 1 || spec.m != 0)
    {
        throw std::invalid_argument(
            name + ": /model/steady_state needs p lags and an intercept and nothing else -- "
                   "p >= 1, n = 1, m = 0; got p = " + std::to_string(spec.p) + ", n = " +
            std::to_string(spec.n) + ", m = " + std::to_string(spec.m));
    }
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword p = static_cast<arma::uword>(spec.p);
    const arma::uword n = k * (k * p + 1);
    if (train.z.n_cols != n)
    {
        throw std::invalid_argument(name + ": /model/steady_state needs /data/train/z of k (k p + 1) "
                                           "= " + std::to_string(n) + " columns, got " +
                                    std::to_string(train.z.n_cols));
    }
    const arma::mat x = compact_regressors(train.z, k);
    if (!arma::all(x.col(x.n_cols - 1) == 1.0))
    {
        throw std::invalid_argument(name + ": /model/steady_state needs the last regressor to be "
                                           "the intercept, a column of ones");
    }
    const arma::uvec lags = steady_state_lag_positions(k, p);
    const arma::uvec intercept = arma::regspace<arma::uvec>(k * k * p, n - 1);
    if (a_prior.v_inv.n_rows == n &&
        arma::abs(arma::mat(a_prior.v_inv.submat(lags, intercept))).max() != 0.0)
    {
        throw std::invalid_argument(
            name + ": /priors/a/v_inv couples the lags with the intercept, and under "
                   "/model/steady_state the intercept has no prior of its own to be coupled with");
    }
    if (mu_prior.mu.n_elem != k || mu_prior.v_inv.n_rows != k || mu_prior.v_inv.n_cols != k)
    {
        throw std::invalid_argument(name + ": /priors/mu needs a mean of k = " + std::to_string(k) +
                                    " and a k x k precision");
    }
    require_finite(mu_prior.mu, "/priors/mu/mu");
    require_finite(mu_prior.v_inv, "/priors/mu/v_inv");
    arma::vec eigenvalues;
    if (!arma::eig_sym(eigenvalues, arma::symmatu(mu_prior.v_inv)) || !(eigenvalues.min() > 0.0))
    {
        throw std::invalid_argument(name + ": /priors/mu/v_inv must be positive definite");
    }
    if (initial_mu.n_elem != k)
    {
        throw std::invalid_argument(name + ": /initial/mu needs k = " + std::to_string(k) +
                                    " elements, got " + std::to_string(initial_mu.n_elem));
    }
    require_finite(initial_mu, "/initial/mu");
}

/// The two Gibbs blocks of the steady-state prior in place of one coefficient
/// draw. `a` and `mu` come in as the last sweep left them and go out drawn;
/// `y` is the stacked response, `z` the SUR regressors (lags then the
/// intercept), `precision` the block-diagonal error precision, (tt k) square,
/// and `a_prior` the base prior, of which the lag block is read. Under
/// VarSpec::stationary the lags are drawn as draw_stationary() draws them, and
/// the return value says whether a stationary draw was found.
inline bool draw_steady_state(arma::vec &a, arma::vec &mu, const arma::vec &y, const arma::mat &z,
                              const arma::sp_mat &precision, const NormalPrior &a_prior,
                              const NormalPrior &mu_prior, const arma::uword k, const arma::uword p,
                              const bool stationary)
{
    const arma::uword tt = y.n_elem / k;
    const arma::uword lag = k * p;
    const arma::mat x = compact_regressors(z, k);
    const arma::mat x_lag = x.cols(0, lag - 1);
    const arma::uvec lags = steady_state_lag_positions(k, p);

    // The lags given mu: the demeaned series on its demeaned lags.
    const arma::rowvec mu_row = arma::repmat(mu.t(), 1, p);
    const arma::mat z_dm = arma::kron(x_lag.each_row() - mu_row, arma::eye<arma::mat>(k, k));
    const arma::vec y_dm = y - arma::repmat(mu, tt, 1);
    const arma::mat v_inv = a_prior.v_inv.submat(lags, lags);
    const arma::mat dz = precision * z_dm;
    const arma::mat post = v_inv + z_dm.t() * dz;
    const arma::vec rhs = v_inv * a_prior.mu.elem(lags) + dz.t() * y_dm;

    arma::vec phi = a.elem(lags);
    bool found = true;
    if (stationary)
    {
        found = draw_stationary(phi, [&]() { return draw_normal_precision(post, rhs); }, k, p);
    }
    else
    {
        phi = draw_normal_precision(post, rhs);
    }

    // mu given the lags: y_t - sum_j A_j y_{t-j} = (I - sum_j A_j) mu + u_t.
    const arma::mat a_lag = arma::reshape(phi, k, lag);
    arma::mat d = arma::eye<arma::mat>(k, k);
    for (arma::uword j = 0; j < p; j++)
    {
        d -= a_lag.cols(j * k, j * k + k - 1);
    }
    const arma::mat w = arma::reshape(y - arma::kron(x_lag, arma::eye<arma::mat>(k, k)) * phi, k, tt);
    arma::mat mu_post = mu_prior.v_inv;
    arma::vec mu_rhs = mu_prior.v_inv * mu_prior.mu;
    for (arma::uword t = 0; t < tt; t++)
    {
        const arma::mat s = arma::mat(precision.submat(t * k, t * k, t * k + k - 1, t * k + k - 1));
        mu_post += d.t() * s * d;
        mu_rhs += d.t() * s * w.col(t);
    }
    mu = draw_normal_precision(mu_post, mu_rhs);

    a.elem(lags) = phi;
    a.tail(k) = d * mu;
    return found;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_STEADY_STATE_SUPPORT_H
