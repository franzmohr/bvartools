// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "kalman_durbin_koopman_2002.h"
#include "bayests/arma.h"

/**
 * @file kalman_durbin_koopman_2002.cpp
 * @brief Simulation smoother of Durbin and Koopman (2002) for a linear Gaussian
 *   state space model.
 */

namespace
{

/// The symmetric square root A of a covariance, so that `A z` carries that
/// covariance for standard normal `z`.
///
/// Deliberately an eigendecomposition and not a Cholesky. Both give an A with
/// `A A' = Sigma`, but they give *different* A, so the same random numbers map
/// to a different draw -- switching would change every posterior this smoother
/// feeds while leaving the distribution intact, which is the most expensive kind
/// of change to verify. A Cholesky would be about three times cheaper per call;
/// that is a trade to make deliberately, not on the way past.
arma::mat symmetric_sqrt(const arma::mat &sigma)
{
  arma::mat U;
  arma::vec s;
  arma::eig_sym(s, U, sigma);
  return U * arma::diagmat(arma::sqrt(s)) * arma::trans(U);
}

} // namespace

/**
 * @brief Draws the state path of a linear Gaussian state space model.
 *
 * Algorithm 2 of Durbin and Koopman (2002) for
 * \f[
 *   y_t = Z_t a_t + u_t, \qquad a_{t+1} = B_t a_t + v_t,
 * \f]
 * with \f$u_t \sim N(0, \Sigma_{u,t})\f$, \f$v_t \sim N(0, \Sigma_{v,t})\f$ and
 * \f$a_0 \sim N(a_{init}, P_{init})\f$, taking into account Jarocinski (2015):
 * the mean of the initial state is set to zero in the first step.
 *
 * @param y K x T matrix of observations, one period per column.
 * @param z KT x M matrix of regressors, the T blocks \f$Z_t\f$ stacked.
 * @param sigma_u the error covariance. Either one K x K matrix for every period,
 *   or a KT x K stack of one per period.
 * @param sigma_v the state innovation covariance. Either one M x M matrix for
 *   every period, or an MT x M stack of one per period.
 * @param B the transition matrix. Either one M x M matrix for every period, or
 *   an MT x M stack of one per period.
 * @param a_init M-vector, the prior mean of \f$a_1\f$ -- the state the *first*
 *   observation loads on, not a period before it. Passed to the filter and only
 *   to the filter: step 1 draws its path from a mean of zero, which is the
 *   Jarocinski (2015) form of the algorithm and is what makes one filtering pass
 *   enough. Both are needed; the mean cannot be carried by step 1 instead, and
 *   `P_init` cannot be dropped -- it is the prior variance of \f$a_1\f$ and
 *   appears in step 1, in the filter's initial P and in the smoothed
 *   \f$a_{init} + P_{init} r_0\f$.
 * @param P_init M x M prior covariance of \f$a_1\f$. Where the caller has
 *   already drawn the previous state, as the TVP samplers do, this is the state
 *   innovation covariance rather than a diffuse prior.
 *
 * @return M x (T+1) matrix of state draws. Column i is the state period i's
 *   observation loads on, for i = 0 ... T-1 -- so the caller wants
 *   `.cols(0, T - 1)`. Column 0 is *not* a state before the sample: it is
 *   \f$a_1\f$, already updated by every observation. Column T is the transition
 *   applied once past the end, informed by no observation, and is there because
 *   the recursions build it on the way; a caller that keeps it has shifted its
 *   whole path by a period.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 */
arma::mat kalman_durbin_koopman_2002(const arma::mat &y, const arma::mat &z,
                                     const arma::mat &sigma_u, const arma::mat &sigma_v,
                                     const arma::mat &B,
                                     const arma::vec &a_init, const arma::mat &P_init) {

  const arma::uword k = y.n_rows;
  const int t = y.n_cols;
  const arma::uword nvars = z.n_cols;

  // Each of the three may arrive as one matrix that holds for every period, or
  // as a stack of one per period -- a time varying error covariance, state
  // innovation covariance and transition are all supported, in any combination.
  //
  // A stride of zero is what lets the rest of the function be written once: the
  // constant form then indexes the same block in every period. The previous
  // spelling instead replicated a constant argument into T copies of itself up
  // front, which allocated a KT x K or MT x M matrix the caller had not asked
  // for and, worse, left the loops below unable to tell that the blocks were
  // identical -- so a constant covariance was eigendecomposed T times.
  const arma::uword u_stride = (sigma_u.n_rows == k) ? 0 : k;
  const arma::uword v_stride = (sigma_v.n_rows == nvars) ? 0 : nvars;
  const arma::uword b_stride = (B.n_rows == nvars) ? 0 : nvars;

  arma::mat yplus = y * 0;
  arma::mat aplus = arma::zeros<arma::mat>(nvars, t + 1);

  // Step 1
  aplus.col(0) = symmetric_sqrt(P_init) * arma::randn<arma::vec>(nvars); // cf. Jarocinski (2015)

  // A covariance that holds for every period has one square root, taken here.
  // A time varying one is decomposed per period inside the loop, as it must be.
  arma::mat A_u, A_v;
  if (u_stride == 0) { A_u = symmetric_sqrt(sigma_u); }
  if (v_stride == 0) { A_v = symmetric_sqrt(sigma_v); }

  for (int i = 0; i < t; i++){
    const arma::uword p1 = i * k;
    const arma::uword p2 = (i + 1) * k - 1;
    const arma::uword up = i * u_stride;
    const arma::uword vp = i * v_stride;
    const arma::uword bp = i * b_stride;

    if (u_stride != 0) { A_u = symmetric_sqrt(sigma_u.rows(up, up + k - 1)); }
    yplus.col(i) = z.rows(p1, p2) * aplus.col(i) + A_u * arma::randn<arma::vec>(k);

    if (v_stride != 0) { A_v = symmetric_sqrt(sigma_v.rows(vp, vp + nvars - 1)); }
    aplus.col(i + 1) = B.rows(bp, bp + nvars - 1) * aplus.col(i) +
                       A_v * arma::randn<arma::vec>(nvars);
  }

  // Kalman filtering
  arma::mat ystar = y - yplus;
  arma::mat a = arma::zeros<arma::mat>(nvars, t + 1);
  a.col(0) = a_init; // cf. DK (2002, p.606)
  arma::mat P = P_init;
  arma::mat v = y * 0;
  arma::mat Fi = arma::zeros<arma::mat>(k * t, k);
  arma::mat K = arma::zeros<arma::mat>(nvars * t, k);
  arma::mat L = arma::zeros<arma::mat>(nvars * t, nvars);
  for (int i = 0; i < t ; i++){
    const arma::uword p1 = i * k;
    const arma::uword p2 = (i + 1) * k - 1;
    const arma::uword pA1 = i * nvars;
    const arma::uword pA2 = (i + 1) * nvars - 1;
    const arma::uword up = i * u_stride;
    const arma::uword vp = i * v_stride;
    const arma::uword bp = i * b_stride;

    const arma::mat B_i = B.rows(bp, bp + nvars - 1);

    v.col(i) = ystar.col(i) - z.rows(p1, p2) * a.col(i);
    Fi.rows(p1, p2) = arma::inv(z.rows(p1, p2) * P * arma::trans(z.rows(p1, p2)) +
                                sigma_u.rows(up, up + k - 1));
    K.rows(pA1, pA2) = B_i * P * arma::trans(z.rows(p1, p2)) * Fi.rows(p1, p2);
    L.rows(pA1, pA2) = B_i - K.rows(pA1, pA2) * z.rows(p1, p2);
    a.col(i + 1) = B_i * a.col(i) + K.rows(pA1, pA2) * v.col(i);
    P = B_i * P * arma::trans(L.rows(pA1, pA2)) + sigma_v.rows(vp, vp + nvars - 1);
  }

  // Backward smoothing
  arma::mat r = arma::zeros<arma::mat>(nvars, t);
  for (int i = (t - 1); i > 0; i--){
    r.col(i - 1) = arma::trans(z.rows(i * k, (i + 1) * k - 1)) * Fi.rows(i * k, (i + 1) * k - 1) * v.col(i) + arma::trans(L.rows(i * nvars, (i + 1) * nvars - 1)) * r.col(i);
  }
  arma::vec r0 = arma::trans(z.rows(0, k - 1)) * Fi.rows(0, k - 1) * v.col(0) + arma::trans(L.rows(0, nvars - 1)) * r.col(0);

  a.col(0) = a_init + P_init * r0;
  for (int i = 0; i < t; i++){
    const arma::uword vp = i * v_stride;
    a.col(i + 1) = B.rows(i * b_stride, i * b_stride + nvars - 1) * a.col(i) +
                   sigma_v.rows(vp, vp + nvars - 1) * r.col(i);
  }

  // Step 3
  return a + aplus;
}
