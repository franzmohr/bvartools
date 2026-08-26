// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef CHAN_JELIAZKOV_2009_H
#define CHAN_JELIAZKOV_2009_H

#include "bayests/arma.h"

arma::mat chan_jeliazkov_2009(const arma::mat &y, const arma::mat &z,
                              const arma::mat &sigma_u, const arma::mat &sigma_v,
                              const arma::mat &B,
                              const arma::vec &a_init, const arma::mat &P_init);

#endif // CHAN_JELIAZKOV_2009_H
