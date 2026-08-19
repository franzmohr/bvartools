// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "core/algorithms/ssvs.h"

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
    blk.order = arma::shuffle(blk.include); // Reorder positions of variable selection

    // Obtain inclusion posterior: the two mixture components evaluated at the
    // current draw, weighted by the prior.
    blk.u0 = 1 / blk.tau0 % arma::exp(-(arma::square(coef) / (2 * blk.tau0sq))) % (1 - blk.inprior);
    blk.u1 = 1 / blk.tau1 % arma::exp(-(arma::square(coef) / (2 * blk.tau1sq))) % blk.inprior;
    blk.post_incl = blk.u1 / (blk.u0 + blk.u1);

    // Draw inclusion parameters in random order
    for (arma::uword i = 0; i < blk.order.n_elem; i++)
    {
        const arma::uword pos = blk.order(i);
        const double draw = (arma::randu() < blk.post_incl(pos)) ? 1.0 : 0.0;
        blk.lambda(pos) = draw;

        // The regressors are untouched; exclusion is expressed as a prior tight
        // enough around zero that the coefficient cannot move.
        prior_v_inv(pos, pos) = 1 / ((draw == 0) ? blk.tau0sq(pos) : blk.tau1sq(pos));
    }
}

} // namespace bayests::core
