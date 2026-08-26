// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_VEC_SUPPORT_H
#define BAYESTS_CORE_MODELS_VEC_SUPPORT_H

#include "bayests/arma.h"

namespace bayests::core
{

/// What separates a VEC from a VAR with wider regressors: the leading k*rank
/// columns of `z` hold `beta' w_{t-1}` kroneckered up, so they are a function of
/// the current draw rather than data, and every VEC sampler here rebuilds them
/// once per iteration. Everything in this file is one of the two readings of
///
///     alpha_t beta_t' w_t,
///
/// linear in alpha or linear in beta, in a constant and a time-varying form.
///
/// `w_t` is k_beta x tt throughout -- the error correction term of one period is
/// a column -- which is the transpose of how TrainData stores it. `beta` is
/// stored as vec of a k_beta x rank matrix and `alpha` as vec of a k x rank one,
/// so the two Kronecker products below are not each other's transpose and the
/// pair is easy to get subtly wrong: `kron(A, B)(i q + a, j p + b) = A(i,j)
/// B(a,b)` is what makes `kron((beta' w)', I_k)` pick out `sum_j (beta' w)_j
/// alpha_ij` and `kron(alpha, w')` pick out `sum_{j,l} alpha_ij w_l beta_lj`,
/// which are the same number.

/// The loadings' regressors when beta does not move: the whole sample in one
/// product, since `beta' w_t` is then r x tt and its transpose kroneckers
/// straight into the (tt k) x (r k) block. The columns past `n_alpha` are data
/// and are left alone.
inline void fill_z_alpha_constant(arma::mat &z, const arma::mat &beta_mat, const arma::mat &w_t,
                                  const int n_alpha, const arma::mat &diag_k)
{
    z.cols(0, n_alpha - 1) = arma::kron(arma::trans(arma::trans(beta_mat) * w_t), diag_k);
}

/// The same when beta is a path: block t of `z` gets kron((beta_t' w_t)', I_k),
/// the k x n_alpha matrix that multiplies vec(alpha_t). `beta` is n_beta x tt.
inline void fill_z_alpha(arma::mat &z, const arma::mat &beta, const arma::mat &w_t, const int k,
                         const int k_beta, const int rank, const arma::mat &diag_k)
{
    const int tt = static_cast<int>(beta.n_cols);
    for (int i = 0; i < tt; i++)
    {
        z.submat(i * k, 0, (i + 1) * k - 1, k * rank - 1) = arma::kron(
            arma::trans(arma::trans(arma::reshape(beta.col(i), k_beta, rank)) * w_t.col(i)),
            diag_k);
    }
}

/// The cointegration vectors' regressors when the loadings do not move: block t
/// of `z_b` gets kron(alpha, w_t'), the k x n_beta matrix that multiplies
/// vec(beta). Period by period even though alpha is constant, because `w_t` is
/// not.
inline void fill_z_beta_constant(arma::mat &z_b, const arma::mat &alpha, const arma::mat &w_t)
{
    const int k = static_cast<int>(alpha.n_rows);
    const int tt = static_cast<int>(w_t.n_cols);
    for (int i = 0; i < tt; i++)
    {
        z_b.rows(i * k, (i + 1) * k - 1) = arma::kron(alpha, arma::trans(w_t.col(i)));
    }
}

/// The same when the loadings are a path: kron(alpha_t, w_t'), with alpha_t read
/// off the front of column t of the n_a x tt coefficient path `a`.
inline void fill_z_beta(arma::mat &z_b, const arma::mat &a, const arma::mat &w_t, const int k,
                        const int rank)
{
    const int tt = static_cast<int>(a.n_cols);
    for (int i = 0; i < tt; i++)
    {
        z_b.rows(i * k, (i + 1) * k - 1) =
            arma::kron(arma::reshape(a.submat(0, i, k * rank - 1, i), k, rank),
                       arma::trans(w_t.col(i)));
    }
}

/// The semi-orthogonal factor of the loadings, alpha (alpha' alpha)^{-1/2}.
///
/// Only the product alpha beta' is identified, so the constant-coefficient VECs
/// split every draw between the two halves: beta is drawn against this
/// normalised Alpha and the scale is handed back afterwards, out of the square
/// root of Beta' Beta. Koop, Leon-Gonzalez and Strachan (2010). The
/// time-varying VECs do not do this -- their state equations carry the
/// normalisation instead, see TvpCointSpacePrior.
inline arma::mat reparameterise_alpha(const arma::mat &alpha, const arma::mat &diag_r)
{
    return alpha * arma::solve(arma::sqrtmat_sympd(arma::trans(alpha) * alpha), diag_r);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_VEC_SUPPORT_H
