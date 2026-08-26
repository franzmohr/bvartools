// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "stochvol_ksc_1998.h"

#include "bayests/arma.h"
#include "core/algorithms/stochvol_mixture.h"

/**
 * @file stochvol_ksc_1998.cpp
 * @brief Posterior draw of a stochastic volatility state with the normal
 *   mixture of Kohn, Shephard and Chib (1998).
 *
 * The mixture is all that distinguishes this from `stochvol_ocsn_2007.cpp`. The
 * draw itself is in `core/algorithms/stochvol_mixture.h`, which says why the two
 * share it rather than each carrying a copy.
 */

namespace
{

/// The seven-component normal mixture of Kohn, Shephard and Chib (1998).
/// Three components fewer than the mixture of Omori, Chib, Shephard and
/// Nakajima (2007) next door, and a coarser approximation to the same log
/// chi-squared distribution for it.
const bayests::core::NormalMixture kMixture = {
    {0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750},
    {-11.40039, -5.24321, -9.83726, 1.50746, -0.65098, 0.52478, -2.35859},
    {5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261}};

} // namespace

/**
 * @brief Draws the log-volatility of a stochastic volatility model.
 *
 * The measurement equation is linearised by squaring and taking logarithms, and
 * the log chi-squared error that leaves is approximated by the seven-component
 * normal mixture of Kohn, Shephard and Chib (1998). Conditional on a mixture
 * indicator per period the state space is linear and Gaussian, so the whole path
 * is drawn in one block from its normal conditional posterior. Each column of
 * `y` is handled independently.
 *
 * See `bayests::core::stochvol_mixture_draw` for the model this belongs to, what
 * each argument has to satisfy, what is thrown when one does not, and the
 * dependence on the global Armadillo random number generator.
 *
 * @param y T x K matrix of error terms, one period per row.
 * @param h T x K matrix of the current log-volatility.
 * @param sigma K-vector of the current variances of the log-volatility
 *   innovations.
 * @param h_init K-vector of the initial states.
 * @param constant K-vector of the offsets added before the logarithm.
 * @return T x K matrix of the drawn log-volatility.
 *
 * @par Example
 * @code
 * arma::mat u = arma::randn<arma::mat>(100, 3);          // T = 100, K = 3
 * arma::mat h = arma::zeros<arma::mat>(100, 3);
 * arma::vec sigma = arma::vec(3).fill(0.05);
 * arma::vec h_init = arma::zeros<arma::vec>(3);
 * arma::vec constant = arma::vec(3).fill(1e-8);
 *
 * h = stochvol_ksc_1998(u, h, sigma, h_init, constant);
 * @endcode
 */
arma::mat stochvol_ksc_1998(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                            const arma::vec &h_init, const arma::vec &constant)
{
  return bayests::core::stochvol_mixture_draw("stochvol_ksc_1998", kMixture, y, h, sigma, h_init,
                                              constant);
}
