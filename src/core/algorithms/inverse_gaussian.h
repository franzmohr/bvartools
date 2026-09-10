// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef INVERSE_GAUSSIAN_H
#define INVERSE_GAUSSIAN_H

#include "bayests/arma.h"

arma::vec inverse_gaussian(const arma::vec &mu, const arma::vec &lambda);

#endif // INVERSE_GAUSSIAN_H
