// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef CONSTRAINED_VAR_PATH_H
#define CONSTRAINED_VAR_PATH_H

#include "bayests/arma.h"
#include "bayests/data.h"

namespace bayests::core
{

/// The Gaussian a VAR puts on its own path, given its parameters: the prior
/// the data-completion step conditions on what was observed.
///
///     y_t = offset_t + A_1 y_{t-1} + ... + A_p y_{t-p} + u_t,   u_t ~ N(0, Sigma_t)
///
/// for t = 0 .. T-1, with the p values before the path known.
///
///   - `presample` is k x p, oldest first: column j is y_{j-p}, so the last
///     column is the period just before the path. Empty when p is zero.
///   - `offset` is k x T, everything in the mean that is not a lag of y --
///     intercepts, trends, exogenous terms -- already multiplied out. Its size
///     is what fixes k and T.
///   - `coefficients` is [A_1 ... A_p], k x kp, the order the lag block of a
///     regressor row has; or T of them stacked, (T k) x kp, one per period.
///   - `covariance` is Sigma_t, k x k; or T of them stacked, (T k) x k.
struct VarPathPrior
{
    arma::mat presample;
    arma::mat offset;
    arma::mat coefficients;
    arma::mat covariance;
};

/// One draw of the path, k x T, from the VAR's prior conditioned on the rows of
/// `constraints` (see bayests::Constraints, with periods counted along the path).
///
/// `soft_variances` holds one error variance per soft group, group g at
/// position g - 1; it may be empty when no row is soft.
///
/// Hard rows of one entry pin that entry, which is then held rather than
/// drawn: the draw is over the entries nothing pins, which is what keeps it
/// cheap when most of a panel is observed. Soft rows enter the precision of
/// that draw. Hard rows of several entries -- an aggregate observed exactly --
/// are imposed afterwards by conditioning the draw on them (Rue and Held 2005,
/// section 2.3.3), so that they hold to rounding in every draw rather than
/// nearly. Everything is done in band storage: the precision is banded because
/// a VAR is Markov of order p and a row reaches over a bounded number of
/// periods, and the band is as wide as the wider of the two.
///
/// Throws std::invalid_argument when the shapes disagree, a position lies
/// outside the path or a soft group has no positive variance, and
/// std::runtime_error when the conditioned system is singular -- hard rows that
/// determine the same combination twice.
arma::mat draw_constrained_var_path(const VarPathPrior &prior, const Constraints &constraints,
                                    const arma::vec &soft_variances);

/// The log density of what the constraints observed, under the VAR's prior,
/// one element per period: element t is the density of the rows whose last
/// period is t given every row that ends before it. Their sum is the log
/// density of all of it, and a period whose rows all end later contributes
/// zero.
///
/// A pinned value counts as observed, so a path pinned whole gives back the
/// VAR's own conditional log likelihood, period by period. That is the
/// quantity a model with a panel not observed whole writes as its pointwise
/// log likelihood: the unobserved entries are integrated out, not filled in.
///
/// Computed by a Kalman filter over a state that carries as many periods as
/// the longer of p and the widest row, which is the prediction error
/// decomposition the elements are.
arma::vec constrained_var_path_log_density(const VarPathPrior &prior,
                                           const Constraints &constraints,
                                           const arma::vec &soft_variances);

/// The mean (k x T) and covariance ((T k) x (T k), in the order of vec of the
/// path) of the distribution draw_constrained_var_path() draws from, from the
/// same factorisation. Dense in the covariance, so for checking the sampler
/// on small paths rather than for a sampler to call.
struct ConstrainedPathMoments
{
    arma::mat mean;
    arma::mat covariance;
};

ConstrainedPathMoments constrained_var_path_moments(const VarPathPrior &prior,
                                                    const Constraints &constraints,
                                                    const arma::vec &soft_variances);

} // namespace bayests::core

#endif // CONSTRAINED_VAR_PATH_H
