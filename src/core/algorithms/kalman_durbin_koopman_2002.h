// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef KALMAN_DURBIN_KOOPMAN_2002_H
#define KALMAN_DURBIN_KOOPMAN_2002_H

#include "bayests/arma.h"

arma::mat kalman_durbin_koopman_2002(const arma::mat &y, const arma::mat &z,
                                     const arma::mat &sigma_u, const arma::mat &sigma_v,
                                     const arma::mat &B,
                                     const arma::vec &a_init, const arma::mat &P_init);

#endif // KALMAN_DURBIN_KOOPMAN_2002_H
