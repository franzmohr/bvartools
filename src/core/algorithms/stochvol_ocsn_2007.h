// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef STOCHVOL_OCSN_2007_H
#define STOCHVOL_OCSN_2007_H

#include "bayests/arma.h"

arma::mat stochvol_ocsn_2007(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                             const arma::vec &h_init, const arma::vec &constant);

#endif // STOCHVOL_OCSN_2007_H
