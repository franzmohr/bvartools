// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_VEC_SUPPORT_H
#define BAYESTS_CORE_MODELS_VEC_SUPPORT_H

#include "bayests/arma.h"
#include "bayests/priors.h"
#include "core/algorithms/truncated_normal.h"

#include <cmath>

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

/// The transition of the cointegration state equation with rho taken out,
/// I_r kron P_tau: the same P_tau for each of the rank relations, which `beta`
/// stacks as vec of a k_beta x rank matrix. An empty `p_tau` is the identity. See
/// TvpCointSpacePrior::p_tau.
inline arma::mat coint_state_transition(const arma::mat &p_tau, const int rank, const int k_beta)
{
    if (p_tau.is_empty())
    {
        const arma::uword n_beta = static_cast<arma::uword>(rank * k_beta);
        return arma::eye<arma::mat>(n_beta, n_beta);
    }
    return arma::kron(arma::eye<arma::mat>(rank, rank), p_tau);
}

/// One draw of rho, the autoregression of the cointegration state equation
///
///     beta_t = rho P beta_{t-1} + eta_t,   eta_t ~ N(0, I),   t = 1, ..., T,
///
/// given the path `beta` (n_beta x tt, one period per column), the state
/// `beta0` of the period before it, and `transition`, the P = I_r kron P_tau of
/// coint_state_transition().
///
/// With the innovation variance fixed at the identity, rho is the coefficient of
/// a regression of beta_t on P beta_{t-1}: the path contributes a normal
/// likelihood in rho whose sufficient statistics are the two sums below,
/// and the prior is uniform on an interval, so the conditional is that normal
/// truncated to the interval -- a Gibbs block, not the Metropolis-within-Gibbs
/// step Koop, Leon-Gonzalez and Strachan (2011) need.
///
/// **The difference is the initial condition, and it is a difference in the
/// model rather than in the algorithm.** Their beta_1 is drawn from the
/// stationary distribution N(0, I / (1 - rho^2)) of the state equation itself,
/// which puts rho in a place no conjugacy survives; here beta_0 has a normal
/// prior of its own, `TvpCointSpacePrior::initial_state`, read from the file and
/// free of rho, so it drops out of this conditional and the draw is exact.
/// Anyone porting their sampler across should know that the prior over the
/// space at the start of the sample is the piece that was not carried over.
///
/// beta_0 does enter through the t = 1 term of the likelihood, which is the
/// whole reason it is an argument.
inline double draw_coint_rho(const arma::mat &beta, const arma::vec &beta0,
                             const arma::mat &transition, const CointRhoPrior &prior)
{
    const arma::uword tt = beta.n_cols;

    // sum_t (P beta_{t-1})' beta_t and sum_t (P beta_{t-1})' (P beta_{t-1}). The
    // lagged path is the path shifted by a column with beta_0 in front, and it is
    // split the way it was before P existed, so that under P = I every sum is
    // formed from the same numbers in the same order.
    const arma::vec lagged0 = transition * beta0;
    double sum_cross = arma::dot(lagged0, beta.col(0));
    double sum_square = arma::dot(lagged0, lagged0);
    if (tt > 1)
    {
        const arma::mat lagged = transition * beta.cols(0, tt - 2);
        sum_cross += arma::accu(lagged % beta.cols(1, tt - 1));
        sum_square += arma::accu(arma::square(lagged));
    }

    // A path that is identically zero says nothing about rho, and its posterior
    // is the prior. Cannot happen to a drawn path, whose innovations have unit
    // variance, but it is what the expressions below divide by.
    if (!(sum_square > 0.0))
    {
        return prior.min + (prior.max - prior.min) * arma::randu<double>();
    }

    return truncated_normal(sum_cross / sum_square, 1.0 / std::sqrt(sum_square), prior.min,
                            prior.max);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_VEC_SUPPORT_H
