// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_VEC_SUPPORT_H
#define BAYESTS_CORE_MODELS_VEC_SUPPORT_H

#include "bayests/arma.h"
#include "bayests/priors.h"
#include "core/algorithms/truncated_normal.h"

#include <cmath>
#include <stdexcept>
#include <string>

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

namespace detail
{

/// The thin SVD x = U diag(s) V' of an n x r matrix with n >= r.
///
/// Every factor this file normalises by is a function of it, and it is taken
/// on x itself rather than by an eigendecomposition of x' x. Squaring x squares
/// its condition number, so a draw that is merely ill-conditioned -- a loading
/// matrix whose singular values span eight orders of magnitude, which a
/// full-rank VEC reaches on real data -- gives an x' x whose smallest
/// eigenvalue is rounding noise. arma::sqrtmat_sympd() refused one that came
/// out below zero and threw partway through a chain; one that came out just
/// above it was inverted into a factor that was not semi-orthogonal at all,
/// without a word. The SVD resolves the same small singular value to working
/// precision.
inline void thin_svd(arma::mat &u, arma::vec &s, arma::mat &v, const arma::mat &x,
                     const char *what)
{
    if (!arma::svd_econ(u, s, v, x))
    {
        throw std::runtime_error(std::string("the singular value decomposition of ") + what +
                                 " failed; the draw is not finite");
    }
}

} // namespace detail

/// The semi-orthogonal factor of the loadings, alpha (alpha' alpha)^{-1/2}.
///
/// Only the product alpha beta' is identified, so the constant-coefficient VECs
/// split every draw between the two halves: beta is drawn against this
/// normalised Alpha and the scale is handed back afterwards by
/// normalise_beta(). Koop, Leon-Gonzalez and Strachan (2010). The
/// time-varying VECs do not do this -- their state equations carry the
/// normalisation instead, see TvpCointSpacePrior.
///
/// With alpha = U S V', alpha (alpha' alpha)^{-1/2} = U S V' V S^{-1} V' = U V',
/// which is what is returned -- see detail::thin_svd() for why it is not
/// computed as written.
inline arma::mat reparameterise_alpha(const arma::mat &alpha)
{
    arma::mat u, v;
    arma::vec s;
    detail::thin_svd(u, s, v, alpha, "the loadings alpha");
    return u * arma::trans(v);
}

/// The other half of that split: the draw `Beta` becomes the semi-orthogonal
/// beta = Beta (Beta' Beta)^{-1/2}, written to `beta`, and its scale
/// (Beta' Beta)^{1/2} is written to `scale` for the caller to hand back to the
/// loadings as Alpha * scale.
///
/// With Beta = U S V', those are U V' and V S V'. The scale is symmetric by
/// construction, which the product of a square root and its inverse was not.
inline void normalise_beta(const arma::mat &Beta, arma::mat &beta, arma::mat &scale)
{
    arma::mat u, v;
    arma::vec s;
    detail::thin_svd(u, s, v, Beta, "the cointegration draw Beta");
    beta = u * arma::trans(v);
    scale = v * arma::diagmat(s) * arma::trans(v);
}

namespace detail
{

/// log |x' P^-1 x| for an n x r matrix x of full column rank, with an empty
/// `p_tau_inv` read as the identity.
///
/// Through the thin SVD x = U S V', as log |S|^2 + log |U' P^-1 U|: the second
/// matrix is as well conditioned as P, so the determinant of the r x r product
/// x' P^-1 x, whose condition number is that of x squared, is never formed. See
/// thin_svd() for why that matters on real data.
inline double log_det_projected(const arma::mat &x, const arma::mat &p_tau_inv, const char *what)
{
    arma::mat u, v;
    arma::vec s;
    thin_svd(u, s, v, x, what);
    double result = 2.0 * arma::accu(arma::log(s));
    if (!p_tau_inv.is_empty())
    {
        double log_det = 0.0;
        if (!arma::log_det_sympd(log_det, arma::mat(arma::symmatu(arma::trans(u) * p_tau_inv * u))))
        {
            throw std::runtime_error(std::string("the cointegration space prior is not positive definite along ") +
                                     what);
        }
        result += log_det;
    }
    return result;
}

} // namespace detail

/// Whether to keep the normal draw `Beta` of the unnormalised cointegration
/// matrix: the Metropolis-Hastings step that makes the constant VECs sample the
/// cointegration space prior they are given when the cointegration term has more
/// rows than the model has equations.
///
/// The prior is Koop, Leon-Gonzalez and Strachan's (2010): beta semi-orthogonal
/// with the matrix angular central Gaussian density |beta' P^-1 beta|^(-k_beta/2),
/// and alpha | beta ~ N(0, v^-1 (beta' P^-1 beta)^-1 kron G). The samplers draw
/// alpha against that, change to A = alpha (alpha' alpha)^(-1/2) and
/// B = beta (alpha' alpha)^(1/2), draw B given A from a normal and split it back
/// with normalise_beta(). Written in A and B, however, the prior is that normal
/// kernel times
///
///     h(B) = |B' P^-1 B|^(-(k_beta - k)/2):
///
/// its own two factors leave |B' P^-1 B|^((k - k_beta)/2) |B' B|^((k_beta - k)/2),
/// and the polar Jacobian from (alpha, beta) to (A, B) contributes
/// |B' B|^((k - k_beta)/2). The factor is one when k_beta = k, the case the paper
/// derives. With deterministic terms restricted to the cointegration space or
/// unmodelled variables in it, k_beta > k, and a normal draw taken as it is
/// overstates |Pi| -- by a quarter to a third of the prior mass in a
/// simulation-based calibration with k = 2 and one or two restricted terms.
///
/// The normal draw is therefore a proposal from the normal part of the
/// conditional, kept with probability min(1, h(Beta) / h(B)) against the current
/// B = beta (alpha' alpha)^(1/2). `alpha` is the current k x rank loading matrix,
/// `beta` the current semi-orthogonal k_beta x rank one. On rejection both stay
/// as they are, which is all the caller has to do. No random number is used when
/// k_beta = k, so the draws of those models do not change.
inline bool accept_coint_draw(const arma::mat &Beta, const arma::mat &alpha, const arma::mat &beta,
                              const arma::mat &p_tau_inv)
{
    const double excess = static_cast<double>(beta.n_rows) - static_cast<double>(alpha.n_rows);
    if (excess <= 0.0)
    {
        return true;
    }

    // |B' P^-1 B| for the current B = beta (alpha' alpha)^(1/2) is
    // |alpha' alpha| |beta' P^-1 beta|, and alpha' alpha is alpha's own Gram
    // matrix, so neither square root is taken.
    const double log_det_proposal = detail::log_det_projected(Beta, p_tau_inv, "the cointegration draw Beta");
    const double log_det_current = detail::log_det_projected(alpha, arma::mat(), "the loadings alpha") +
                                   detail::log_det_projected(beta, p_tau_inv, "the cointegration matrix beta");
    const double log_ratio = -0.5 * excess * (log_det_proposal - log_det_current);
    if (std::isnan(log_ratio))
    {
        throw std::runtime_error("the acceptance ratio of the cointegration draw is not a number");
    }
    return log_ratio >= 0.0 || std::log(arma::randu<double>()) < log_ratio;
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
