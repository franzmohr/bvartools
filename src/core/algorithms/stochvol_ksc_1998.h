// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef STOCHVOL_KSC_1998_H
#define STOCHVOL_KSC_1998_H

#include "bayests/arma.h"

arma::mat stochvol_ksc_1998(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                             const arma::vec &h_init, const arma::vec &constant);

#endif // STOCHVOL_KSC_1998_H
