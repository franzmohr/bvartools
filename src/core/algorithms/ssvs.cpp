// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "core/algorithms/ssvs.h"

#include "core/algorithms/inclusion_probability.h"

#include <cmath>

namespace bayests::core
{

SsvsBlock::SsvsBlock(const arma::vec &initial_lambda, const VarSelPrior &prior)
    : lambda(initial_lambda),
      inprior(prior.inprior),
      tau0(prior.ssvs.tau0),
      tau1(prior.ssvs.tau1),
      tau0sq(arma::square(prior.ssvs.tau0)),
      tau1sq(arma::square(prior.ssvs.tau1)),
      include(prior.include)
{
}

void ssvs_sweep(SsvsBlock &blk, const arma::vec &coef, arma::mat &prior_v_inv)
{
    for (arma::uword i = 0; i < blk.include.n_elem; i++)
    {
        const arma::uword pos = blk.include(i);
        const double square = coef(pos) * coef(pos);

        // The two mixture components evaluated at the current draw and weighted
        // by the prior, in logs. Formed as densities, both underflow to zero for
        // a coefficient many slab widths from zero, and their ratio is then a
        // NaN that excludes the one coefficient the data most want in.
        const double l1 = std::log(blk.inprior(pos)) - std::log(blk.tau1(pos)) -
                          square / (2 * blk.tau1sq(pos));
        const double l0 = std::log(1 - blk.inprior(pos)) - std::log(blk.tau0(pos)) -
                          square / (2 * blk.tau0sq(pos));

        const bool included = arma::randu() < inclusion_probability(l1 - l0);
        blk.lambda(pos) = included ? 1.0 : 0.0;

        // The regressors are untouched; exclusion is expressed as a prior tight
        // enough around zero that the coefficient cannot move.
        prior_v_inv(pos, pos) = 1 / (included ? blk.tau1sq(pos) : blk.tau0sq(pos));
    }
}

} // namespace bayests::core
