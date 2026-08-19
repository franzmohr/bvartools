// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_BVS_H
#define BAYESTS_CORE_ALGORITHMS_BVS_H

#include "bayests/arma.h"
#include "bayests/priors.h"

#include <cmath>

namespace bayests::core
{

/// Bayesian variable selection, as a Gibbs block.
///
/// One inclusion indicator per selected coefficient, swept in a fresh random
/// order every iteration. The indicator masks the regressors themselves and is
/// accepted or rejected on a likelihood ratio, so exclusion is exact -- a
/// coefficient that is out contributes nothing. That is the difference from
/// SSVS, in ssvs.h, which leaves the regressors alone and instead moves the
/// prior precision between a spike and a slab; the two are alternatives, chosen
/// by VarSelection, and no sampler runs both at once.
///
/// A sampler holds one block per coefficient vector it selects over -- in
/// practice one for the lag coefficients `a` and one for the contemporaneous
/// coefficients `psi`. That is what this type is for: the state a sweep carries
/// between draws used to be a dozen loose locals per block, distinguished only
/// by an `a_` or `psi_` prefix.
///
/// Every RNG call the scheme makes happens inside bvs_sweep(), in the order the
/// samplers consumed them before this was factored out. That ordering is
/// load-bearing: the regression harness in test/ pins the seed and fingerprints
/// the output, so a draw added, dropped or reordered here shows up as a changed
/// fingerprint on every fixture that exercises selection.

/// What a position stands for when BVS switches it off.
enum class BvsScope
{
    /// One coefficient, addressed linearly. This is the whole story where the
    /// coefficients are a vector, one value per regressor.
    element,

    /// One regressor across the whole sample: the coefficients are a path,
    /// nparams x tt, and switching a regressor off zeroes its row in every
    /// period.
    path_row,
};

/// State a BVS sweep carries between draws.
struct BvsBlock
{
    /// Inclusion indicators as they stood at the end of the last sweep. Read by
    /// the prescreen and only refreshed once the sweep is over, so a decision
    /// taken early in a sweep does not change the prescreen for a position
    /// visited later in the same sweep.
    arma::vec lambda;

    /// diagmat(lambda), which is what actually masks the regressors and the
    /// coefficients. This is the copy a sweep updates as it goes.
    arma::mat lambda_diag;

    arma::vec lprior_0; ///< log(1 - prior inclusion probability)
    arma::vec lprior_1; ///< log(prior inclusion probability)

    /// Positions the sweep may touch, and those positions reshuffled -- kept as
    /// a member only so the permutation does not reallocate every draw.
    arma::uvec include;
    arma::uvec order;

    BvsBlock() = default;
    BvsBlock(const arma::vec &initial_lambda, const VarSelPrior &prior);
};

namespace detail
{

/// Switching a position off and back on, for both shapes of coefficients.
///
/// `.row()` works on a vector as well as a matrix, and for a vector the two
/// scopes coincide -- a one-column matrix has its linear index equal to its row
/// index -- so a single pair of templates covers every caller.
template <class Coefficients>
inline void switch_position_off(Coefficients &theta, const arma::uword pos, const BvsScope scope)
{
    if (scope == BvsScope::path_row)
    {
        theta.row(pos).zeros();
    }
    else
    {
        theta(pos) = 0.0;
    }
}

template <class Coefficients>
inline void switch_position_on(Coefficients &theta, const Coefficients &coef,
                               const arma::uword pos, const BvsScope scope)
{
    if (scope == BvsScope::path_row)
    {
        theta.row(pos) = coef.row(pos);
    }
    else
    {
        theta(pos) = coef(pos);
    }
}

} // namespace detail

/// One BVS sweep: visit the selected positions in random order and accept or
/// reject each one on a likelihood ratio against the current mask.
///
/// `coef` is masked in place on the way out, so it holds the drawn coefficients
/// with the excluded positions zeroed -- a vector for a constant-coefficient
/// model, an nparams x tt path for a time-varying one.
///
/// `log_likelihood(candidate)` must return the log likelihood of a candidate
/// coefficient set up to the prior, that is -r'(Omega)r / 2 for whatever
/// residual `candidate` implies and whatever precision the block is conditioned
/// on. Keeping it a callable is what lets one sweep serve every model: the
/// residual is a single matrix product for a constant-coefficient model and a
/// loop over periods for a time-varying one, and the precision is dense in some
/// models and sparse in others. The sweep adds the log prior itself.
///
/// Two candidates are evaluated per surviving position -- the mask with this
/// position forced off, and with it forced on -- against the mask as it stood
/// before the sweep began.
template <class Coefficients, class LogLikelihood>
void bvs_sweep(BvsBlock &blk, Coefficients &coef, const BvsScope scope,
               LogLikelihood &&log_likelihood)
{
    blk.order = arma::shuffle(blk.include);

    const Coefficients masked = blk.lambda_diag * coef;
    Coefficients candidate;

    for (arma::uword j = 0; j < blk.order.n_elem; j++)
    {
        const arma::uword pos = blk.order(j);

        // Prescreen: a position only gets the full likelihood treatment with
        // probability equal to the prior of the state it is currently in. This
        // is the one draw that happens whether or not the position is
        // revisited, so it is taken before anything else.
        const double lprior_current =
            (blk.lambda(pos) == 1) ? blk.lprior_1(pos) : blk.lprior_0(pos);
        if (std::log(arma::randu()) >= lprior_current)
        {
            continue;
        }

        candidate = masked;
        detail::switch_position_off(candidate, pos, scope);
        const double l0 = log_likelihood(candidate) + blk.lprior_0(pos);

        candidate = masked;
        detail::switch_position_on(candidate, coef, pos, scope);
        const double l1 = log_likelihood(candidate) + blk.lprior_1(pos);

        blk.lambda_diag(pos, pos) = (l1 - l0 >= std::log(arma::randu())) ? 1.0 : 0.0;
    }

    coef = blk.lambda_diag * coef;
    blk.lambda = arma::vectorise(blk.lambda_diag.diag());
}

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_BVS_H
