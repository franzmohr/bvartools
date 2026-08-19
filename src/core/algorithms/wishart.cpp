// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "wishart.h"
#include "bayests/arma.h"

/**
 * @file wishart.cpp
 * @brief Posterior draw of a precision matrix from its Wishart conditional.
 */

/**
 * @brief Draws a K x K precision matrix from its Wishart posterior.
 *
 * Given a series of K-dimensional error terms \f$u_t\f$, \f$t = 1, \dots, T\f$,
 * with
 * \f[
 *   u_t \sim N(0, \Sigma),
 * \f]
 * and an inverse Wishart prior \f$\Sigma \sim IW(S_0, \nu_0)\f$ on the
 * covariance matrix, the conditional posterior of the precision matrix
 * \f$\Sigma^{-1}\f$ is Wishart with scale matrix
 * \f[
 *   \left( S_0 + \sum_{t=1}^{T} u_t u_t' \right)^{-1}
 *      = \left( S_0 + U U' \right)^{-1}
 * \f]
 * and \f$\nu_0 + T\f$ degrees of freedom, where \f$U = [u_1, \dots, u_T]\f$.
 * The function draws one sample from that distribution.
 *
 * Note that the degrees of freedom are *not* derived from the data inside the
 * function. The caller must pass the already updated posterior value
 * `post_df`, which is usually the prior degrees of freedom plus the number of
 * observations \f$T\f$.
 *
 * @param u K x T matrix of error terms. Each column contains the K errors of
 *   one period, so that `u * u.t()` is the K x K matrix of sums of squared
 *   errors. Taken by const reference to avoid a copy, so a caller may pass a
 *   temporary such as `arma::reshape(y - z * a, k, tt)` directly. The caller's
 *   matrix is free to change between calls; it is only the draw that treats it
 *   as fixed.
 * @param prior_scale K x K symmetric positive definite prior scale matrix
 *   \f$S_0\f$ of the inverse Wishart prior on the covariance matrix. Taken by
 *   const reference for the same reason.
 * @param post_df Posterior degrees of freedom \f$\nu_0 + T\f$. Must be at
 *   least K for the draw to be positive definite with probability one.
 *
 * @return K x K draw of the precision matrix \f$\Sigma^{-1}\f$. Invert it, for
 *   example with `arma::inv_sympd`, to obtain the corresponding covariance
 *   matrix.
 *
 * @pre `prior_scale + u * u.t()` is symmetric positive definite. The
 *   inversion uses `arma::inv_sympd` and only reads the lower triangle, so a
 *   non-symmetric argument is silently symmetrised rather than rejected.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 *
 * @par Example
 * @code
 * arma::mat u = arma::randn<arma::mat>(3, 100);   // K = 3, T = 100
 * arma::mat prior_scale = arma::eye<arma::mat>(3, 3);
 * int post_df = 3 + 100;                          // prior df + T
 *
 * arma::mat sigma_inv = wishart(u, prior_scale, post_df);
 * arma::mat sigma = arma::inv_sympd(sigma_inv);
 * @endcode
 *
 * @see arma::wishrnd
 */
arma::mat wishart(const arma::mat &u, const arma::mat &prior_scale, int post_df) {
  return arma::wishrnd(arma::inv_sympd(prior_scale + u * u.t()), post_df);
}
