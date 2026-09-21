// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef STOCHVOL_OCSN_2007_H
#define STOCHVOL_OCSN_2007_H

#include "bayests/arma.h"

arma::mat stochvol_ocsn_2007(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                             const arma::vec &h_init, const arma::vec &constant);

/// The non-centred draw with the same mixture: see
/// `bayests::core::stochvol_mixture_draw_noncentred`.
arma::mat stochvol_ocsn_2007_noncentred(const arma::mat &y, const arma::mat &h,
                                        const arma::vec &h_init, const arma::vec &omega,
                                        const arma::vec &constant, arma::mat &y_centred,
                                        arma::mat &precision);

#endif // STOCHVOL_OCSN_2007_H
