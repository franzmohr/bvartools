// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "core/algorithms/bvs.h"

namespace bayests::core
{

BvsBlock::BvsBlock(const arma::vec &initial_lambda, const VarSelPrior &prior)
    : lambda(initial_lambda),
      lambda_diag(arma::diagmat(initial_lambda)),
      lprior_0(arma::log(1 - prior.inprior)),
      lprior_1(arma::log(prior.inprior)),
      include(prior.include)
{
}

} // namespace bayests::core
