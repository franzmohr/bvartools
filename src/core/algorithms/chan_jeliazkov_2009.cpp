// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "chan_jeliazkov_2009.h"

#include "bayests/arma.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @file chan_jeliazkov_2009.cpp
 * @brief Precision based draw of a state path, after Chan and Jeliazkov (2009).
 *
 * The alternative to `kalman_durbin_koopman_2002` for the same conditional
 * posterior. Rather than filtering forward and sampling backward, the whole path
 * is treated as one Gaussian vector whose precision is block tridiagonal --
 * because the state equation is first order Markov and the measurement of one
 * period touches one state -- and drawn in a single pass over a block banded
 * Cholesky. The two are the same distribution reached two ways, so which is
 * faster is a question about constants rather than about orders.
 *
 * Both take their arguments in the same shapes and return the same thing, so a
 * caller can swap one for the other, and a test can put the same inputs through
 * both and compare what comes out.
 */

namespace
{

void require(bool ok, const std::string &what)
{
  if (!ok) {
    throw std::invalid_argument("chan_jeliazkov_2009: " + what);
  }
}

std::string dims(const arma::mat &m)
{
  return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

/// True if every element is a finite number, read off the exponent bits.
///
/// Not `arma::is_finite()` and not `std::isfinite()`, for the reason
/// `stochvol_mixture.h` sets out at length: a host that compiles these sources
/// with `-ffast-math` is licensed to fold either of them to `true`.
bool all_finite(const arma::mat &m)
{
  static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

  const arma::uword n = m.n_elem;
  for (arma::uword i = 0; i < n; i++) {
    std::uint64_t bits;
    const double value = m[i];
    std::memcpy(&bits, &value, sizeof(bits));
    if ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL) {
      return false;
    }
  }
  return true;
}

/// The inverse of a covariance, by the symmetric positive definite route.
arma::mat precision_of(const arma::mat &sigma, const char *what)
{
  arma::mat out;
  require(arma::inv_sympd(out, sigma),
          std::string(what) + " is not symmetric positive definite, so it has no precision");
  return out;
}

} // namespace

/**
 * @brief Draws the state path of a linear Gaussian state space model.
 *
 * For
 * \f[
 *   y_t = Z_t s_t + u_t, \qquad s_t = \sum_{j=1}^{p} A_{j,t} s_{t-j} + v_t,
 * \f]
 * with \f$u_t \sim N(0, \Sigma_{u,t})\f$, \f$v_t \sim N(0, \Sigma_{v,t})\f$ and
 * the first p states drawn jointly from \f$N(a_{init}, P_{init})\f$, the whole
 * path \f$s = (s_0', \dots, s_T')'\f$ is Gaussian, and its precision K is block
 * banded of bandwidth p. Three things populate it, and nothing else does:
 *
 * - the prior, contributing \f$P_{init}^{-1}\f$ across the first p block rows
 *   and columns, and \f$P_{init}^{-1}a_{init}\f$ to the first p blocks of b;
 * - the measurement of period t, contributing
 *   \f$Z_t'\Sigma_{u,t}^{-1}Z_t\f$ to block (t, t) and
 *   \f$Z_t'\Sigma_{u,t}^{-1}y_t\f$ to block t of b -- one state per observation,
 *   so this touches the diagonal only;
 * - the transition producing state t, whose residual
 *   \f$s_t - \sum_j A_{j,t} s_{t-j}\f$ contributes \f$C_a' \Sigma_{v,t}^{-1}
 *   C_b\f$ to block (t-a, t-b) for a, b = 0 ... p, with \f$C_0 = I\f$ and
 *   \f$C_j = -A_{j,t}\f$. Those indices span t-p to t, which is where the
 *   bandwidth comes from: a first order transition couples neighbours, an order
 *   p one couples p apart, and no arrangement of the arguments couples further.
 *
 * The mean is \f$K^{-1}b\f$. The Cholesky factor of a block banded matrix is
 * block banded of the same width, so it comes out of one sweep over the periods
 * whatever p is, and the T(T+1)M^2 matrix is never formed.
 *
 * Indexing follows `kalman_durbin_koopman_2002`: observation `i` is measured
 * against state column `i`, and the transition producing column `t` is the one
 * indexed `t - 1`. The last column has no observation of its own.
 *
 * @param y K x T matrix of observations, one period per column.
 * @param z KT x M matrix of regressors, the T blocks \f$Z_t\f$ stacked.
 * @param sigma_u the error covariance. Either one K x K matrix for every period,
 *   or a KT x K stack of one per period.
 * @param sigma_v the state innovation covariance. Either one M x M matrix for
 *   every period, or an MT x M stack of one per period.
 * @param B the transition, of order p. The p coefficient matrices side by side,
 *   \f$[A_1 \; A_2 \; \dots \; A_p]\f$, so M x pM for a transition that holds
 *   for every period or MT x pM for a stack of one per period. p is read off the
 *   width, so an M x M argument is the first order case and everything below
 *   reduces to it. A stacked argument's row block `t - 1` is the one that
 *   produces column `t`; its leading `p - 1` row blocks are never read, since no
 *   transition produces a column the prior already covers.
 * @param a_init pM-vector, the prior mean of the first p states, in time order:
 *   \f$(a_1', a_2', \dots, a_p')'\f$. For p = 1 this is the state the first
 *   observation loads on.
 * @param P_init pM x pM prior covariance of those first p states, jointly. Must
 *   be positive definite; see the note on that below.
 *
 * @return M x (T+1) matrix of state draws. Column i is the state period i's
 *   observation loads on, for i = 0 ... T-1 -- so the caller wants
 *   `.cols(0, T - 1)`. Column 0 is *not* a state before the sample: it is
 *   \f$a_1\f$, conditioned on every observation. Column T is the transition
 *   applied once past the end, informed by no observation.
 *
 * @note `P_init` has to be invertible here, where `kalman_durbin_koopman_2002`
 *   also accepts a singular one. A precision based sampler needs
 *   \f$P_{init}^{-1}\f$, and the limiting case -- a state fixed at `a_init`
 *   rather than merely tightly distributed around it -- is a different model,
 *   one state block shorter. Use the simulation smoother for it, or a small
 *   `P_init`.
 *
 * @throws std::invalid_argument if the dimensions of the arguments do not
 *   describe a state space model, if any of them is not finite, or if `P_init`
 *   or a covariance has no symmetric positive definite inverse. The message
 *   names the argument.
 * @throws std::runtime_error if the posterior precision turns out not to be
 *   positive definite, or if the drawn path is not finite.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 */
arma::mat chan_jeliazkov_2009(const arma::mat &y, const arma::mat &z,
                              const arma::mat &sigma_u, const arma::mat &sigma_v,
                              const arma::mat &B,
                              const arma::vec &a_init, const arma::mat &P_init)
{
  const arma::uword k = y.n_rows;
  const arma::uword tt = y.n_cols;
  const arma::uword m = z.n_cols;

  // Checks
  require(k > 0 && tt > 0, "'y' must have at least one row and one column, got " + dims(y));
  require(m > 0, "'z' must have at least one column");
  require(z.n_rows == k * tt,
          "'z' must have as many rows as 'y' has elements (" + std::to_string(k * tt) +
              "), got " + std::to_string(z.n_rows));
  // The order of the transition is read off the width of B: p coefficient
  // matrices side by side. One of them is the first order case.
  require(B.n_cols > 0 && B.n_cols % m == 0,
          "'B' must have a multiple of " + std::to_string(m) + " columns, one block per lag, got " +
              std::to_string(B.n_cols));
  const arma::uword p = B.n_cols / m;
  require(tt + 1 >= p, "a transition of order " + std::to_string(p) + " needs at least that many "
                       "state columns, and " + std::to_string(tt) + " periods give " +
                       std::to_string(tt + 1));

  require(a_init.n_elem == p * m, "'a_init' must have " + std::to_string(p * m) +
                                      " elements, one per state the prior covers, got " +
                                      std::to_string(a_init.n_elem));
  require(P_init.n_rows == p * m && P_init.n_cols == p * m,
          "'P_init' must be " + std::to_string(p * m) + "x" + std::to_string(p * m) + ", got " +
              dims(P_init));
  require(all_finite(y), "'y' contains NaN or infinite values");
  require(all_finite(z), "'z' contains NaN or infinite values");
  require(all_finite(a_init), "'a_init' contains NaN or infinite values");

  // Each of the three may arrive as one matrix that holds for every period, or
  // as a stack of one per period. A stride of zero is what lets the body be
  // written once, and it is also what lets a constant argument be inverted once
  // rather than T times -- the same arrangement `kalman_durbin_koopman_2002`
  // uses, for the same reason.
  const arma::uword u_stride = (sigma_u.n_rows == k) ? 0 : k;
  const arma::uword v_stride = (sigma_v.n_rows == m) ? 0 : m;
  const arma::uword b_stride = (B.n_rows == m) ? 0 : m;

  require(u_stride == 0 || sigma_u.n_rows == k * tt,
          "'sigma_u' must be " + std::to_string(k) + "x" + std::to_string(k) + " or " +
              std::to_string(k * tt) + "x" + std::to_string(k) + ", got " + dims(sigma_u));
  require(v_stride == 0 || sigma_v.n_rows == m * tt,
          "'sigma_v' must be " + std::to_string(m) + "x" + std::to_string(m) + " or " +
              std::to_string(m * tt) + "x" + std::to_string(m) + ", got " + dims(sigma_v));
  require(b_stride == 0 || B.n_rows == m * tt,
          "'B' must have " + std::to_string(m) + " or " + std::to_string(m * tt) +
              " rows, got " + std::to_string(B.n_rows));
  require(sigma_u.n_cols == k, "'sigma_u' must have " + std::to_string(k) + " columns, got " +
                                   std::to_string(sigma_u.n_cols));
  require(sigma_v.n_cols == m, "'sigma_v' must have " + std::to_string(m) + " columns, got " +
                                   std::to_string(sigma_v.n_cols));

  // A constant covariance is inverted here, once. A time varying one is
  // inverted per period below, as it must be.
  arma::mat u_prec, v_prec;
  if (u_stride == 0) { u_prec = precision_of(sigma_u, "'sigma_u'"); }
  if (v_stride == 0) { v_prec = precision_of(sigma_v, "'sigma_v'"); }

  const arma::mat p_init_prec = precision_of(P_init, "'P_init'");

  // Every time varying parameter model in this library follows a random walk, so
  // B is the identity and both B'Sigma_v^-1 B and B'Sigma_v^-1 are Sigma_v^-1.
  // Tested once here, because otherwise they are two M x M x M products per
  // period spent multiplying by one -- which at M = 45 is as much arithmetic as
  // the factorisation itself.
  const bool b_is_identity =
      p == 1 && b_stride == 0 &&
      arma::approx_equal(B, arma::eye<arma::mat>(m, m), "absdiff", 0.0);

  // The band, upper triangle only, since K is symmetric: `block(i, d)` is the
  // M x M block at row i and column i + d, for d = 0 ... p. Overwritten by the
  // Cholesky factor in place, which has the same shape, so the whole draw holds
  // (T+1)(p+1) matrices of M x M and the T(T+1)M^2 matrix never exists.
  const arma::uword n = tt + 1;
  std::vector<arma::mat> band(n * (p + 1), arma::zeros<arma::mat>(m, m));
  const auto block = [&band, p](arma::uword i, arma::uword d) -> arma::mat & {
    return band[i * (p + 1) + d];
  };
  arma::mat rhs(m, n, arma::fill::zeros);

  // The prior, over the first p states jointly. For p = 1 this is one block.
  for (arma::uword i = 0; i < p; i++) {
    for (arma::uword j = i; j < p; j++) {
      block(i, j - i) +=
          p_init_prec.submat(i * m, j * m, (i + 1) * m - 1, (j + 1) * m - 1);
    }
    rhs.col(i) += p_init_prec.rows(i * m, (i + 1) * m - 1) * a_init;
  }

  // The measurements. Period t loads on state column t alone, so this reaches
  // the main diagonal only, whatever p is.
  for (arma::uword t = 0; t < tt; t++) {
    if (u_stride != 0) {
      u_prec = precision_of(sigma_u.rows(t * u_stride, t * u_stride + k - 1), "'sigma_u'");
    }
    const arma::mat Z_t = z.rows(t * k, (t + 1) * k - 1);
    const arma::mat zu = Z_t.t() * u_prec;
    block(t, 0) += zu * Z_t;
    rhs.col(t) += zu * y.col(t);
  }

  // The transitions. The one producing state t is indexed t - 1, and exists for
  // t = p ... T -- the first p states are the prior's, not any transition's.
  //
  // Writing C_0 = I and C_j = -A_j, the residual is sum_a C_a s_{t-a} and its
  // quadratic form puts C_a' V C_b in block (t-a, t-b). Only a >= b is stored,
  // that being the upper triangle, and the rest is its transpose.
  for (arma::uword t = p; t <= tt; t++) {
    const arma::uword bi = (t - 1) * b_stride;
    const arma::uword vi = (t - 1) * v_stride;
    if (v_stride != 0) {
      v_prec = precision_of(sigma_v.rows(vi, vi + m - 1), "'sigma_v'");
    }

    block(t, 0) += v_prec; // a = b = 0

    if (b_is_identity) {
      // The random walk, spelled out: A_1 = I makes both products below the
      // identity's, and skipping them is worth about a tenth of the assembly.
      block(t - 1, 1) += -v_prec;
      block(t - 1, 0) += v_prec;
      continue;
    }

    const arma::mat B_t = B.rows(bi, bi + m - 1);
    for (arma::uword j = 1; j <= p; j++) {
      const arma::mat aj = B_t.cols((j - 1) * m, j * m - 1);
      const arma::mat ajv = aj.t() * v_prec;

      // a = j, b = 0: block (t - j, t).
      block(t - j, j) += -ajv;

      // a = j, b = l, keeping a >= b so that t - a <= t - b.
      for (arma::uword l = 1; l <= j; l++) {
        block(t - j, j - l) += ajv * B_t.cols((l - 1) * m, l * m - 1);
      }
    }
  }

  // Block Cholesky, K = R'R with R block upper bidiagonal: R_i on the diagonal
  // and S_i above it, from
  //
  //   R_0'R_0 = D_0,   R_i'S_i = U_i,   R_{i+1}'R_{i+1} = D_{i+1} - S_i'S_i,
  //
  // one sweep, each step a dense M x M factorisation or triangular solve. This
  // is the whole reason the algorithm is worth having: the same posterior as a
  // T(T+1)M^2 dense factorisation, at O(T M^3).
  //
  // K is positive definite by construction -- a sum of a positive definite
  // prior precision and positive semidefinite quadratic forms -- so a failure
  // means a conditioning value has already degenerated.
  //
  // Written out as it stands rather than in terms of the structure the blocks
  // often have, because the structured spelling was tried and measured slower.
  // With a random walk over a diagonal state covariance -- every time varying
  // parameter model here -- U_i is diagonal and the Schur complement needs no S:
  //
  //   S_i'S_i = U_i' (R_i'R_i)^-1 U_i = Lam_i D_i^-1 Lam_i,
  //
  // a symmetric inverse scaled by a vector either side, which on paper replaces
  // a triangular solve with M right hand sides and a full M x M x M product by
  // M^2 work. On this machine it came out 15 to 30% *slower* at every size
  // measured, because `inv_sympd` factorises the block a second time and the
  // inverse it then computes is a poorer LAPACK kernel than the product it
  // replaced -- fewer flops at a worse rate. The dense spelling below hands the
  // largest term to `gemm`, which is the fastest kernel available, and that wins.
  // From (R'R)_{i,j} = sum_k R_{k,i}' R_{k,j}, where R_{k,i} is non-zero only for
  // i - p <= k <= i, the k that contribute to block (i, j) run from j - p to i.
  // Splitting off k = i,
  //
  //   R_{i,i}' R_{i,j} = K_{i,j} - sum_{k=j-p}^{i-1} R_{k,i}' R_{k,j},
  //
  // which is a Cholesky at j = i and a triangular solve above it. One sweep, and
  // for p = 1 it is the bidiagonal recursion it generalises.
  // The accumulation is done into the band itself rather than into a scratch
  // matrix: block (i, d) holds K_{i,i+d} until step i consumes it, and is never
  // read again afterwards, so subtracting in place costs nothing and saves an
  // M x M copy per block -- which at M = 45 was a quarter of this loop.
  for (arma::uword i = 0; i < n; i++) {
    for (arma::uword d = 0; d <= p && i + d < n; d++) {
      const arma::uword j = i + d;
      for (arma::uword q = (j >= p) ? j - p : 0; q < i; q++) {
        block(i, d) -= block(q, i - q).t() * block(q, j - q);
      }

      arma::mat r;
      if (d == 0) {
        if (!arma::chol(r, block(i, 0))) {
          throw std::runtime_error("chan_jeliazkov_2009: the posterior precision of the state path "
                                   "is not positive definite at period " + std::to_string(i));
        }
      } else if (!arma::solve(r, arma::trimatl(block(i, 0).t()), block(i, d),
                              arma::solve_opts::fast)) {
        throw std::runtime_error("chan_jeliazkov_2009: the Cholesky factor of the state path "
                                 "precision is singular at period " + std::to_string(i));
      }
      block(i, d) = r;
    }
  }

  // Forward substitution, R'f = b, and then the backward pass. The mean and the
  // noise share it: R x = f + w solves both at once, because R m = f and
  // R n = w make R(m + n) = f + w, and Cov(R^-1 w) = (R'R)^-1 = K^-1.
  arma::mat f(m, n);
  for (arma::uword i = 0; i < n; i++) {
    arma::vec b = rhs.col(i);
    for (arma::uword q = (i >= p) ? i - p : 0; q < i; q++) {
      b -= block(q, i - q).t() * f.col(q);
    }
    arma::vec solved;
    if (!arma::solve(solved, arma::trimatl(block(i, 0).t()), b, arma::solve_opts::fast)) {
      throw std::runtime_error("chan_jeliazkov_2009: the state path posterior is singular at "
                               "period " + std::to_string(i));
    }
    f.col(i) = solved;
  }

  arma::mat draw(m, n);
  for (arma::uword i = n; i-- > 0;) {
    arma::vec b = f.col(i) + arma::randn<arma::vec>(m);
    for (arma::uword d = 1; d <= p && i + d < n; d++) {
      b -= block(i, d) * draw.col(i + d);
    }
    arma::vec solved;
    if (!arma::solve(solved, arma::trimatu(block(i, 0)), b, arma::solve_opts::fast)) {
      throw std::runtime_error("chan_jeliazkov_2009: the state path posterior is singular at "
                               "period " + std::to_string(i));
    }
    draw.col(i) = solved;
  }

  if (!all_finite(draw)) {
    throw std::runtime_error("chan_jeliazkov_2009: the drawn state path is not finite");
  }

  return draw;
}
