// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef WISHART_H
#define WISHART_H

#include "bayests/arma.h"

arma::mat wishart(const arma::mat &u, const arma::mat &prior_scale, int post_df);

#endif // WISHART_H
