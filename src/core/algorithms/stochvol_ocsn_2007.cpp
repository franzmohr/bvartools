// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "stochvol_ocsn_2007.h"

#include "bayests/arma.h"
#include "core/algorithms/stochvol_mixture.h"

/**
 * @file stochvol_ocsn_2007.cpp
 * @brief Posterior draw of a stochastic volatility state with the normal
 *   mixture of Omori, Chib, Shephard and Nakajima (2007).
 *
 * The mixture is all that distinguishes this from `stochvol_ksc_1998.cpp`. The
 * draw itself is in `core/algorithms/stochvol_mixture.h`, which says why the two
 * share it rather than each carrying a copy.
 */

namespace
{

/// The ten-component normal mixture of Omori, Chib, Shephard and Nakajima
/// (2007). Three components more than the mixture of Kim, Shephard and Chib
/// (1998) next door, and a closer approximation to the same log chi-squared
/// distribution for it.
const bayests::core::NormalMixture kMixture = {
    {0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115},
    {1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384,
     -14.65000},
    {0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342}};

} // namespace

/**
 * @brief Draws the log-volatility of a stochastic volatility model.
 *
 * The measurement equation is linearised by squaring and taking logarithms, and
 * the log chi-squared error that leaves is approximated by the ten-component
 * normal mixture of Omori, Chib, Shephard and Nakajima (2007). Conditional on a
 * mixture indicator per period the state space is linear and Gaussian, so the
 * whole path is drawn in one block from its normal conditional posterior. Each
 * column of `y` is handled independently.
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
 * h = stochvol_ocsn_2007(u, h, sigma, h_init, constant);
 * @endcode
 */
arma::mat stochvol_ocsn_2007(const arma::mat &y, const arma::mat &h, const arma::vec &sigma,
                             const arma::vec &h_init, const arma::vec &constant)
{
  return bayests::core::stochvol_mixture_draw("stochvol_ocsn_2007", kMixture, y, h, sigma, h_init,
                                              constant);
}
