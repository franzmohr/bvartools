// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef KALMAN_DURBIN_KOOPMAN_2002_H
#define KALMAN_DURBIN_KOOPMAN_2002_H

#include "bayests/arma.h"

arma::mat kalman_durbin_koopman_2002(arma::mat &y, arma::mat &z,
                                     arma::mat sigma_u, arma::mat sigma_v,
                                     arma::mat B,
                                     arma::vec &a_init, arma::mat &P_init);

#endif // KALMAN_DURBIN_KOOPMAN_2002_H
