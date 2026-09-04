// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef CHAN_JELIAZKOV_2009_H
#define CHAN_JELIAZKOV_2009_H

#include "bayests/arma.h"

arma::mat chan_jeliazkov_2009(const arma::mat &y, const arma::mat &z,
                              const arma::mat &sigma_u, const arma::mat &sigma_v,
                              const arma::mat &B,
                              const arma::vec &a_init, const arma::mat &P_init);

/// The same draw with the trailing `known.n_rows` elements of every state column
/// held at the values in `known` rather than drawn. Returns one column per
/// period, and no trailing column to drop. See the source for what that means
/// for `z` and `y`, which are unchanged from the whole-state form.
arma::mat chan_jeliazkov_2009_conditional(const arma::mat &y, const arma::mat &z,
                                          const arma::mat &sigma_u, const arma::mat &sigma_v,
                                          const arma::mat &B,
                                          const arma::vec &a_init, const arma::mat &P_init,
                                          const arma::mat &known);

#endif // CHAN_JELIAZKOV_2009_H
