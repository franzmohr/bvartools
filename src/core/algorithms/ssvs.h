// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_SSVS_H
#define BAYESTS_CORE_ALGORITHMS_SSVS_H

#include "bayests/arma.h"
#include "bayests/priors.h"

namespace bayests::core
{

/// Stochastic search variable selection, as a Gibbs block.
///
/// One inclusion indicator per selected coefficient, swept in a fresh random
/// order every iteration. The regressors are left alone; the indicator moves the
/// prior precision between a spike and a slab, so exclusion is not exact but a
/// prior tight enough around zero that the coefficient cannot move. That is the
/// difference from BVS, in bvs.h, which masks the regressors and so excludes
/// exactly; the two are alternatives, chosen by VarSelection, and no sampler
/// runs both at once.
///
/// Offered only by the two constant-coefficient models with a normal prior on
/// the coefficients -- VarNormalGamma and VarNormalWishart. The stochastic
/// volatility and time-varying parameter samplers reject it in validate(),
/// because a spike-and-slab prior on a coefficient that is itself a random walk
/// would have to be a prior on the state equation rather than on the level.
///
/// Every RNG call the scheme makes happens inside ssvs_sweep(), in the order the
/// samplers consumed them before this was factored out. That ordering is
/// load-bearing: the regression harness in test/ pins the seed and fingerprints
/// the output, so a draw added, dropped or reordered here shows up as a changed
/// fingerprint on every fixture that exercises selection.

/// State an SSVS sweep carries between draws.
struct SsvsBlock
{
    arma::vec lambda;
    arma::vec inprior;

    arma::vec tau0; ///< Spike standard deviation
    arma::vec tau1; ///< Slab standard deviation
    arma::vec tau0sq;
    arma::vec tau1sq;

    /// Positions the sweep may touch, and those positions reshuffled -- kept as
    /// a member only so the permutation does not reallocate every draw.
    arma::uvec include;
    arma::uvec order;

    /// Mixture weights and the resulting inclusion posterior. Members rather
    /// than locals for the same reason as `order`: they are the same size every
    /// draw.
    arma::vec u0;
    arma::vec u1;
    arma::vec post_incl;

    SsvsBlock() = default;
    SsvsBlock(const arma::vec &initial_lambda, const VarSelPrior &prior);
};

/// One SSVS sweep: draw every inclusion indicator from its posterior and point
/// the prior precision at the matching mixture component.
///
/// `prior_v_inv` is the prior precision of `coef` and is modified in place --
/// only the diagonal entries at the selected positions, which is why the
/// samplers keep a mutable copy of the prior rather than referring to the input.
void ssvs_sweep(SsvsBlock &blk, const arma::vec &coef, arma::mat &prior_v_inv);

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_SSVS_H
