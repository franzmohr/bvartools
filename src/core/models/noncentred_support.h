// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_NONCENTRED_SUPPORT_H
#define BAYESTS_CORE_MODELS_NONCENTRED_SUPPORT_H

#include "bayests/priors.h"
#include "bayests/results.h"
#include "core/algorithms/kalman_durbin_koopman_2002.h"
#include "core/algorithms/stochvol_ocsn_2007.h"
#include "core/models/model_support.h"

#include <cmath>
#include <stdexcept>
#include <string>

/**
 * @file noncentred_support.h
 * @brief The non-centred random walk of Frühwirth-Schnatter and Wagner (2010),
 *   and the Savage-Dickey ordinates Chan (2018) tests time variation with.
 *
 * A random walk with a start and an innovation variance,
 * \f[
 *   x_t = x_{t-1} + v_t, \qquad v_t \sim N(0, \sigma), \qquad x_0 \text{ given},
 * \f]
 * is the same model as
 * \f[
 *   x_t = x_0 + \omega \tilde x_t, \qquad
 *   \tilde x_t = \tilde x_{t-1} + e_t, \qquad e_t \sim N(0, 1), \qquad \tilde x_0 = 0,
 * \f]
 * with \f$\omega = \pm\sqrt{\sigma}\f$. Written this way, \f$\omega\f$ enters the
 * measurement equation as a regression coefficient beside \f$x_0\f$, so a normal
 * prior \f$\omega \sim N(0, V_\omega)\f$ is conjugate, and \f$\omega = 0\f$ -- the
 * state does not move -- is an interior point of its support instead of the
 * boundary of the variance's. That is what the Savage-Dickey density ratio
 * needs: the Bayes factor of the time-varying model against the one with the
 * state held at \f$x_0\f$ is
 * \f[
 *   BF = \frac{p(\omega = 0)}{p(\omega = 0 \mid y)},
 * \f]
 * and the denominator is estimated by averaging, over the draws, the density at
 * zero of \f$\omega\f$'s conditional posterior given the standardised path. That
 * conditional is normal, so the ordinate is exact and each draw contributes one
 * number. The numerator is the prior's and needs no draws. The implied prior
 * on \f$\sigma = \omega^2\f$ is Gamma(1/2, 1 / (2 V_\omega)), which puts more mass
 * near zero than the usual inverse gamma does.
 *
 * The sign of \f$\omega\f$ is not identified -- flipping it together with the
 * whole of \f$\tilde x\f$ leaves the likelihood where it was -- so its posterior
 * is symmetric, and bimodal wherever the data want the state to move. Each draw
 * flips the pair with probability one half, as Frühwirth-Schnatter and Wagner
 * suggest, so the chain visits both modes rather than one. The path
 * \f$x_t\f$, and \f$\sigma = \omega^2\f$ with it, is untouched by the flip, and so
 * is the ordinate at zero.
 *
 * What is written out is the log of the ordinate, per state and for the whole
 * block jointly. The per-state ones test one state at a time. The joint one
 * tests "every state moves" against "none does" -- not "some state moves" -- so
 * a block in which one state moves and many stand still can come out against
 * time variation, each state that stands still costing about a log point. It
 * is not the sum of the per-state ones either, since the regression couples
 * the states' \f$\omega\f$ through the data. Chan (2018, section 2.3) warns that
 * the estimate is noisy exactly where the Bayes factor is large -- the ordinate
 * is then an average of tiny numbers -- so a host should report a numerical
 * standard error beside it.
 */

namespace bayests::core
{

/// One draw of \f$(x_0, \omega)\f$ from its joint normal conditional posterior,
/// and the log ordinates of \f$\omega\f$'s marginal at zero.
struct NoncentredCoefficients
{
    arma::vec x0;    ///< n; the state before the sample.
    arma::vec omega; ///< n; the signed standard deviations.

    /// n; \f$\log p(\omega_i = 0 \mid \cdot)\f$, \f$x_0\f$ and the other
    /// \f$\omega_j\f$ integrated out.
    arma::vec log_zero;

    /// \f$\log p(\omega = 0 \mid \cdot)\f$ for the whole block, \f$x_0\f$
    /// integrated out.
    double log_zero_joint = 0.0;
};

/// Draws \f$(x_0, \omega)\f$ given what the data say about them and the prior.
///
/// Conditional on the standardised path, the measurement equation is a linear
/// regression on the 2n coefficients \f$(x_0', \omega')'\f$, so the caller sums
/// its Gram matrix and cross product -- `data_precision` and `data_rhs`, in that
/// order of coefficients -- and this adds the two independent priors, \f$x_0 \sim
/// N(\mu, V^{-1})\f$ and \f$\omega_i \sim N(0, V_{\omega,i})\f$.
///
/// Independence of the two priors is also the condition the Savage-Dickey
/// ratio needs: the prior of \f$x_0\f$ under the constant model has to be the
/// conditional one under the time-varying model at \f$\omega = 0\f$.
///
/// @param data_precision 2n x 2n; \f$\sum_t X_t' S_t X_t\f$.
/// @param data_rhs 2n; \f$\sum_t X_t' S_t y_t\f$.
/// @param x0_prior prior on the state before the sample.
/// @param omega_v n; prior variances of \f$\omega\f$.
inline NoncentredCoefficients draw_noncentred_coefficients(const arma::mat &data_precision,
                                                           const arma::vec &data_rhs,
                                                           const NormalPrior &x0_prior,
                                                           const arma::vec &omega_v)
{
    const arma::uword n = omega_v.n_elem;
    const double log_2pi = std::log(2.0 * arma::datum::pi);

    arma::mat precision = data_precision;
    precision.submat(0, 0, n - 1, n - 1) += x0_prior.v_inv;
    for (arma::uword i = 0; i < n; i++)
    {
        precision(n + i, n + i) += 1.0 / omega_v(i);
    }
    precision = arma::symmatu(precision);

    arma::vec rhs = data_rhs;
    rhs.head(n) += x0_prior.v_inv * x0_prior.mu;

    // The ordinates need omega's marginal covariance, which is a block of the
    // inverse rather than of the precision. 2n is the number of states, not of
    // periods, so the inverse is cheap next to the path draw that precedes it.
    arma::mat covariance;
    if (!arma::inv_sympd(covariance, precision))
    {
        throw std::runtime_error("the posterior precision of the non-centred state equation is "
                                 "not positive definite");
    }
    const arma::vec mean = covariance * rhs;
    const arma::vec omega_mean = mean.tail(n);
    const arma::mat omega_cov = arma::symmatu(covariance.submat(n, n, 2 * n - 1, 2 * n - 1));

    NoncentredCoefficients out;

    const arma::vec omega_var = omega_cov.diag();
    out.log_zero = -0.5 * (log_2pi + arma::log(omega_var) + arma::square(omega_mean) / omega_var);

    double log_det = 0.0;
    if (!arma::log_det_sympd(log_det, omega_cov))
    {
        throw std::runtime_error("the posterior covariance of the non-centred standard deviations "
                                 "is not positive definite");
    }
    out.log_zero_joint = -0.5 * (static_cast<double>(n) * log_2pi + log_det +
                                 arma::dot(omega_mean, arma::solve(omega_cov, omega_mean)));

    const arma::vec draw = draw_normal_precision(precision, rhs);
    out.x0 = draw.head(n);
    out.omega = draw.tail(n);

    return out;
}

/// Flips each \f$(\omega_i, \tilde x_i)\f$ pair with probability one half.
///
/// `x_tilde` holds state i in row i when `states_in_rows`, and in column i
/// otherwise -- the coefficient paths are laid out one way and the
/// log-volatility the other.
inline void switch_noncentred_signs(arma::vec &omega, arma::mat &x_tilde, const bool states_in_rows)
{
    const arma::vec u = arma::randu<arma::vec>(omega.n_elem);
    for (arma::uword i = 0; i < omega.n_elem; i++)
    {
        if (u(i) < 0.5)
        {
            omega(i) = -omega(i);
            if (states_in_rows)
            {
                x_tilde.row(i) *= -1.0;
            }
            else
            {
                x_tilde.col(i) *= -1.0;
            }
        }
    }
}

/// Stores one draw's \f$\omega\f$ and ordinates at column `draw_pos`.
inline void store_noncentred(NoncentredStateDraws &out, const arma::uword draw_pos,
                             const arma::vec &omega, const NoncentredCoefficients &c)
{
    out.omega.col(draw_pos) = omega;
    out.log_zero.col(draw_pos) = c.log_zero;
    out.log_zero_joint(0, draw_pos) = c.log_zero_joint;
}

/// Sizes the three matrices of NoncentredStateDraws for n states.
inline void allocate_noncentred(NoncentredStateDraws &out, const arma::uword n,
                                const arma::uword iterations)
{
    out.omega = arma::mat(n, iterations);
    out.log_zero = arma::mat(n, iterations);
    out.log_zero_joint = arma::mat(1, iterations);
}

/// The whole non-centred step for a block whose states are regression
/// coefficients: the coefficients of a VAR, or the covariance block.
///
/// Draws the standardised path given \f$(x_0, \omega)\f$ with the simulation
/// smoother, then \f$(x_0, \omega)\f$ given the path, then the signs, and
/// rebuilds the path \f$x_t = x_0 + \omega \tilde x_t\f$ the rest of the sampler
/// reads. Measurement and state equation are
/// \f[
///   y_t - Z_t x_0 = Z_t \operatorname{diag}(\omega)\, \tilde x_t + u_t, \qquad
///   \tilde x_t = \tilde x_{t-1} + e_t, \qquad \tilde x_0 = 0,
/// \f]
/// so the smoother runs with a unit state variance and \f$\tilde x_1 \sim N(0,
/// I)\f$, and the regression for \f$(x_0, \omega)\f$ has the per-period design
/// \f$X_t = [Z_t, Z_t \operatorname{diag}(\tilde x_t)]\f$.
///
/// @param y r x T responses, one period per column.
/// @param z rT x n regressors, the T blocks \f$Z_t\f$ stacked, already masked
///   by any selection the caller performs.
/// @param covariance_blocks rT x r error covariances, one block per period.
/// @param precision_blocks rT x r their inverses.
/// @param prior the block's prior; `omega_v` must be set.
/// @param x0 n; updated.
/// @param omega n; updated.
/// @param x_tilde n x T standardised path; overwritten.
/// @param path n x T; overwritten with \f$x_0 + \omega \tilde x_t\f$.
/// @return the draw of \f$(x_0, \omega)\f$ with its ordinates at zero.
inline NoncentredCoefficients draw_noncentred_path(const arma::mat &y, const arma::mat &z,
                                                   const arma::mat &covariance_blocks,
                                                   const arma::mat &precision_blocks,
                                                   const RandomWalkPrior &prior, arma::vec &x0,
                                                   arma::vec &omega, arma::mat &x_tilde,
                                                   arma::mat &path)
{
    const arma::uword r = y.n_rows;
    const arma::uword tt = y.n_cols;
    const arma::uword n = z.n_cols;

    // The standardised path, given where the state starts and how far it moves.
    arma::mat y_offset(r, tt);
    for (arma::uword t = 0; t < tt; t++)
    {
        y_offset.col(t) = y.col(t) - z.rows(t * r, (t + 1) * r - 1) * x0;
    }
    arma::mat z_scaled = z;
    z_scaled.each_row() %= arma::trans(omega);
    const arma::mat identity = arma::eye<arma::mat>(n, n);
    x_tilde = kalman_durbin_koopman_2002(y_offset, z_scaled, covariance_blocks, identity, identity,
                                         arma::zeros<arma::vec>(n), identity)
                  .cols(0, tt - 1);

    // Where it starts and how far it moves, given the standardised path.
    arma::mat data_precision = arma::zeros<arma::mat>(2 * n, 2 * n);
    arma::vec data_rhs = arma::zeros<arma::vec>(2 * n);
    arma::mat x_t(r, 2 * n);
    for (arma::uword t = 0; t < tt; t++)
    {
        const arma::mat z_t = z.rows(t * r, (t + 1) * r - 1);
        x_t.cols(0, n - 1) = z_t;
        x_t.cols(n, 2 * n - 1) = z_t.each_row() % arma::trans(x_tilde.col(t));
        const arma::mat sx = precision_blocks.rows(t * r, (t + 1) * r - 1) * x_t;
        data_precision += arma::trans(x_t) * sx;
        data_rhs += arma::trans(sx) * y.col(t);
    }

    NoncentredCoefficients c =
        draw_noncentred_coefficients(data_precision, data_rhs, prior.initial_state, prior.omega_v);
    x0 = c.x0;
    omega = c.omega;

    switch_noncentred_signs(omega, x_tilde, true);

    path = x_tilde.each_col() % omega;
    path.each_col() += x0;

    return c;
}

/// The whole non-centred step for the log-volatility of a stochastic volatility
/// model, one state per variable.
///
/// After the mixture approximation the measurement equation of variable i is
/// \f$\log(u_{it}^2 + c_i) - m_{s_t} = h_{i0} + \omega_i \tilde h_{it} +
/// \varepsilon_{it}\f$ with \f$\varepsilon_{it} \sim N(0, v_{s_t})\f$, so given the
/// indicators and \f$\tilde h\f$ it is a regression on \f$(h_{i0}, \omega_i)\f$,
/// and the k regressions are stacked into one because the prior on the initial
/// log-volatilities may couple them. The data precision is block diagonal
/// across variables, so the per-variable ordinates are also the conditionally
/// independent ones and the joint ordinate is their sum wherever that prior is
/// diagonal.
///
/// @param u T x k orthogonalised errors.
/// @param prior the log-volatility's state prior; `omega_v` must be set.
/// @param offset k; added inside the logarithm.
/// @param h T x k log-volatility; the conditioning value of the indicators on
///   the way in, the new draw on the way out.
/// @param h_init k; updated.
/// @param omega k; updated.
/// @param h_tilde T x k standardised log-volatility; overwritten.
inline NoncentredCoefficients draw_noncentred_log_volatility(const arma::mat &u,
                                                             const RandomWalkPrior &prior,
                                                             const arma::vec &offset, arma::mat &h,
                                                             arma::vec &h_init, arma::vec &omega,
                                                             arma::mat &h_tilde)
{
    const arma::uword k = u.n_cols;

    arma::mat y_centred, precision;
    h_tilde = stochvol_ocsn_2007_noncentred(u, h, h_init, omega, offset, y_centred, precision);

    arma::mat data_precision = arma::zeros<arma::mat>(2 * k, 2 * k);
    arma::vec data_rhs = arma::zeros<arma::vec>(2 * k);
    for (arma::uword i = 0; i < k; i++)
    {
        const arma::vec p = precision.col(i);
        const arma::vec p_h = p % h_tilde.col(i);
        data_precision(i, i) = arma::accu(p);
        data_precision(i, k + i) = arma::accu(p_h);
        data_precision(k + i, i) = data_precision(i, k + i);
        data_precision(k + i, k + i) = arma::dot(p_h, h_tilde.col(i));
        data_rhs(i) = arma::dot(p, y_centred.col(i));
        data_rhs(k + i) = arma::dot(p_h, y_centred.col(i));
    }

    NoncentredCoefficients c =
        draw_noncentred_coefficients(data_precision, data_rhs, prior.initial_state, prior.omega_v);
    h_init = c.x0;
    omega = c.omega;

    switch_noncentred_signs(omega, h_tilde, false);

    h = h_tilde.each_row() % arma::trans(omega);
    h.each_row() += arma::trans(h_init);

    return c;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_NONCENTRED_SUPPORT_H
