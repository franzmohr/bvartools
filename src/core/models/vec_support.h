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

/// What the draw of the cointegration matrix conditions on: see
/// augment_loadings().
struct CointDrawLoadings
{
    /// k x rank: the loadings' rows of the semi-orthogonal A. The data enter the
    /// draw of B through them, and the loadings after the draw are these rows
    /// times the scale normalise_beta() returns. With k_beta = k they are the
    /// Alpha of reparameterise_alpha().
    arma::mat top;
    /// n_beta x n_beta prior precision of vec(B), the k_beta x rank draw.
    arma::mat prior_vinv;
};

/// The collapsed Gibbs step of Koop, Leon-Gonzalez and Strachan (2010) for a
/// cointegration term with any number of rows: the loadings' side of the draw of
/// the cointegration matrix, and the exact normal prior that draw is taken from.
///
/// The prior is theirs: beta semi-orthogonal with the matrix angular central
/// Gaussian density |beta' P^-1 beta|^(-k_beta/2), and
/// alpha | beta ~ N(0, v^-1 (beta' P^-1 beta)^-1 kron G). The samplers draw
/// alpha against it, change to A = alpha (alpha' alpha)^(-1/2) and
/// B = beta (alpha' alpha)^(1/2), draw B given A from a normal and split it back
/// with normalise_beta(). That B is normal given A is the paper's Proposition 1,
/// which needs alpha and beta to have the same number of rows. With
/// deterministic terms restricted to the cointegration space or unmodelled
/// variables in it, k_beta > k, and in (A, B) the prior is then the normal kernel
/// times |B' P^-1 B|^(-(k_beta - k)/2): the exponents of the density of beta and
/// of the normaliser of alpha | beta no longer cancel, nor do the two polar
/// Jacobians. A normal draw taken as it is overstates |Pi|.
///
/// Rather than correct the normal draw -- as a Metropolis-Hastings proposal it
/// was kept so rarely on the country models of a global VEC, where k_beta - k is
/// four to seven, that a chain could stay at its starting values for thousands
/// of draws -- the loadings are given the rows they lack. k_beta - k auxiliary
/// rows are drawn from
///
///     alpha_aux | beta ~ N(0, c^-1 (beta' Q beta)^-1 kron gamma^-1 I),
///
/// independent of alpha, the other coefficients, the error precision and the
/// data, so the model for everything else is unchanged. Stacked beneath alpha
/// they make a k_beta x rank matrix alpha_+, and the change to
/// A = alpha_+ (alpha_+' alpha_+)^(-1/2) and B = beta (alpha_+' alpha_+)^(1/2)
/// is the paper's with equal dimensions:
///
/// - With v > 0, c = v and Q = P^-1, so alpha_+ | beta is the paper's prior with
///   G_+ = diag(G, gamma^-1 I). The normaliser of alpha_+ | beta cancels the
///   density of beta, the polar Jacobians cancel each other, and the prior in
///   (A, B) is the normal kernel with precision
///   v (A' G_+^-1 A) kron P^-1.
/// - With v = 0 the prior on alpha is flat and the one on the space uniform,
///   whatever P is, which is how the samplers have always read a zero
///   shrinkage. Then c = 1 and Q = I: the normaliser of alpha_aux is one for a
///   semi-orthogonal beta, and the prior in (A, B) is the normal kernel with
///   precision gamma (A_aux' A_aux) kron I.
///
/// Pi = alpha beta' = A_top B', with A_top the first k rows of A, so the data
/// enter the draw of B through A_top and the posterior of B given A is normal
/// exactly -- no correction, no rejection. After the draw, alpha = A_top times
/// the scale of B, and the auxiliary rows are discarded: they are drawn afresh
/// from their conditional before every draw of B, which is what lets the other
/// blocks go on sampling the model without them.
///
/// `g_inv` is the precision G^-1 the loadings' prior uses. gamma is its mean
/// diagonal element, which puts the auxiliary rows on the scale of the loadings.
/// Any positive gamma leaves the posterior unchanged and could only affect how
/// fast the chain mixes; on the US model of test/unit_coint_dees_us.cpp, scaling
/// it by anything from 1e-4 to 1e4 did not measurably. `p_tau_inv` may be empty, read as the identity.
///
/// With k_beta = k there is nothing to augment: no random number is drawn, and
/// `top` and `prior_vinv` are exactly what the samplers used before, so the draws
/// of those models do not change.
inline CointDrawLoadings augment_loadings(const arma::mat &alpha, const arma::mat &beta,
                                          const arma::mat &g_inv, const double v_inv,
                                          const arma::mat &p_tau_inv)
{
    const arma::uword k = alpha.n_rows;
    const arma::uword rank = alpha.n_cols;
    const arma::uword k_beta = beta.n_rows;
    const arma::mat p_inv =
        p_tau_inv.is_empty() ? arma::mat(arma::eye<arma::mat>(k_beta, k_beta)) : p_tau_inv;

    CointDrawLoadings out;

    if (k_beta <= k)
    {
        out.top = reparameterise_alpha(alpha);
        out.prior_vinv = arma::kron(arma::trans(out.top) * g_inv * out.top, v_inv * p_inv);
        return out;
    }

    const arma::uword n_aux = k_beta - k;
    const double gamma = arma::trace(g_inv) / static_cast<double>(k);
    if (!std::isfinite(gamma) || gamma <= 0.0)
    {
        throw std::runtime_error("the error precision is not positive, so the auxiliary loadings of "
                                 "the cointegration draw have no scale");
    }

    const arma::mat space =
        v_inv > 0.0 ? arma::mat(v_inv * p_inv) : arma::mat(arma::eye<arma::mat>(k_beta, k_beta));

    // alpha_aux = Z C^-T / sqrt(gamma) with Z standard normal and C' C the row
    // precision c beta' Q beta: vec(Z M) has covariance M' M kron I, and
    // M' M = C^-1 C^-T = (C' C)^-1.
    arma::mat chol_upper;
    if (!arma::chol(chol_upper, arma::mat(arma::symmatu(arma::trans(beta) * space * beta))))
    {
        throw std::runtime_error("the cointegration space prior is not positive definite along the "
                                 "cointegration matrix beta");
    }
    const arma::mat aux =
        arma::trans(arma::solve(arma::trimatu(chol_upper),
                                arma::mat(arma::trans(arma::randn<arma::mat>(n_aux, rank))))) /
        std::sqrt(gamma);

    const arma::mat a_full = reparameterise_alpha(arma::join_cols(alpha, aux));
    out.top = a_full.rows(0, k - 1);
    const arma::mat a_aux = a_full.rows(k, k_beta - 1);

    out.prior_vinv = arma::kron(arma::trans(out.top) * g_inv * out.top, v_inv * p_inv) +
                     arma::kron(gamma * (arma::trans(a_aux) * a_aux), space);
    return out;
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
