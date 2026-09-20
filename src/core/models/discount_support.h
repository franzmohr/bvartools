// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_DISCOUNT_SUPPORT_H
#define BAYESTS_CORE_MODELS_DISCOUNT_SUPPORT_H

#include "bayests/arma.h"
#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/reporter.h"
#include "core/models/model_support.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace bayests::core
{

/// The recursion behind every discounted model here: the matrix normal dynamic
/// linear model of West and Harrison (1997, ch. 16) with Uhlig's (1997)
/// discounted Wishart on the error precision.
///
/// It is a header of its own because a VAR and a VEC with a fixed cointegration
/// space differ in the design matrix they hand it and in nothing else. The
/// design is `tt` rows by `n_design` columns whatever the model calls those
/// columns, every equation shares them -- which is the whole reason the
/// posterior is conjugate -- and what a column means is the caller's business.
///
/// Both models' in-sample score and their score past the end of the sample run
/// the same DiscountState::step(), so the two cannot drift apart into
/// describing different models from the same file.

/// Log density of the k dimensional Student t with `nu` degrees of freedom and
/// scale matrix `psi`, at the forecast error `e`.
///
/// Returns a quiet NaN rather than throwing when the scale is not positive
/// definite: one unusable period should not lose a whole path that is otherwise
/// fine, and the caller reports it as a missing score.
inline double log_mvt(const arma::vec &e, const arma::mat &psi, const double nu)
{
    const double k = static_cast<double>(e.n_elem);

    arma::mat chol_psi;
    if (!arma::chol(chol_psi, arma::symmatu(psi), "lower"))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double log_det = 2.0 * arma::accu(arma::log(chol_psi.diag()));
    const arma::vec sol = arma::solve(arma::trimatl(chol_psi), e);

    return std::lgamma(0.5 * (nu + k)) - std::lgamma(0.5 * nu) -
           0.5 * k * std::log(nu * arma::datum::pi) - 0.5 * log_det -
           0.5 * (nu + k) * std::log1p(arma::dot(sol, sol) / nu);
}

/// The marginal scale of every coefficient. Theta_t is matrix t about m with
/// left covariance C and right covariance S, so element (i, j) has scale
/// sqrt(C_ii S_jj). Returned as vec(B_t) with B_t of k x n_design, which is
/// row-major in Theta and the ordering every other model in this project
/// vectorises a coefficient matrix by.
inline arma::vec coefficient_scale(const arma::mat &c, const arma::mat &s)
{
    const arma::vec c_diag = c.diag();
    const arma::vec s_diag = s.diag();
    return arma::vectorise(arma::sqrt(c_diag * s_diag.t()).t());
}

/// One draw of Sigma, scaled so that drawing the coefficients against it
/// reproduces the matrix t the filter reports.
///
/// `t_df` is the degrees of freedom the filter carries, which are the Student t
/// ones: the predictive is T_{t_df} and `coefficient_scale()` reports
/// sqrt(C_ii S_jj) as the scale of element (i, j). The inverse Wishart needs
/// different ones. With Sigma ~ IW_k(Psi, nu) and vec(Theta) | Sigma normal
/// with covariance Sigma kron C, element (i, j) of Theta is marginally t on
/// **nu - k + 1** degrees of freedom with scale squared
/// C_ii Psi_jj / (nu - k + 1). The filter's point estimate is S = Psi / t_df,
/// so Psi = t_df S, and requiring the t degrees of freedom to come out at
/// `t_df` gives nu = t_df + k - 1, whereupon the scale squared is exactly
/// C_ii S_jj.
///
/// Passing `t_df` straight through as nu instead is the mistake this converts
/// away: it leaves the draws a factor sqrt(t_df / (t_df - k + 1)) too wide,
/// which is 6.6 per cent at three variables and a discount of 0.94 -- small
/// enough to read as Monte Carlo error and uniform across every coefficient,
/// which is what gives it away.
///
/// `chol_inv_scale` is the lower Cholesky factor of (t_df S)^-1, the inverse of
/// Psi.
inline arma::mat draw_inverse_wishart(const arma::mat &chol_inv_scale, const double t_df,
                                      const arma::uword k)
{
    const double nu = t_df + static_cast<double>(k) - 1.0;

    // Bartlett: L A with A lower triangular, chi on the diagonal and standard
    // normals below it, gives a Wishart draw's factor.
    arma::mat a(k, k, arma::fill::zeros);
    for (arma::uword i = 0; i < k; i++)
    {
        a(i, i) = std::sqrt(arma::chi2rnd(nu - static_cast<double>(i)));
        for (arma::uword j = 0; j < i; j++)
        {
            a(i, j) = arma::randn();
        }
    }
    const arma::mat factor = chol_inv_scale * a;
    return arma::inv_sympd(arma::symmatu(factor * factor.t()));
}

/// The response as one row per period, in whichever of the three layouts
/// TrainData allows it to arrive in.
///
/// Through stacked_response() rather than by vectorising `train.y` directly,
/// because a single row or a single column holding vec(y') is how the HDF5
/// files store it and only the stacked vector is the same object in all three
/// cases. Vectorising a period-per-row matrix gives vec(y), not vec(y'), and so
/// transposes the sample without a word.
inline arma::mat response_by_period(const TrainData &train, const arma::uword k,
                                    const arma::uword tt)
{
    return arma::reshape(stacked_response(train), k, tt).t();
}

/// The posterior as the filter carries it: (Theta | Sigma) matrix normal about
/// `m` with left covariance `c`, and Sigma inverse Wishart with `df` degrees of
/// freedom and sum of squares `d`, whose point estimate is `d / df`.
///
/// A plain struct rather than a class with a constructor, because it is started
/// from a prior in one place and from a stored posterior in another: a score
/// past the end of the sample resumes the same filter at the last period the
/// sample pinned down.
struct DiscountState
{
    arma::mat m;     ///< n_design x k.
    arma::mat c;     ///< n_design x n_design.
    double df = 0.0; ///< Degrees of freedom of the inverse Wishart.
    arma::mat d;     ///< k x k sum of squares; the scale is `d / df`.

    /// Evolves by the two discounts, scores `y` against the prediction, and
    /// updates on it. Returns the log predictive density of the period and
    /// writes the predictive mean to `mean` where one is wanted.
    ///
    /// The evolution is a model rather than a plug-in: `R = C / delta` is
    /// exactly the predicted covariance of a random walk whose innovation
    /// covariance is `((1 - delta) / delta) C`, and discounting the Wishart
    /// leaves its point estimate alone while letting it drift.
    double step(const arma::vec &z, const arma::vec &y, const double d_beta,
                const double d_sigma, arma::vec *mean = nullptr)
    {
        const arma::mat r = c / d_beta;
        const double df_pred = d_sigma * df;
        const arma::mat d_pred = d_sigma * d;

        // Because the equations share the design, the whole measurement
        // covariance is q Sigma with q a scalar: no k x k solve here, and no
        // (n_design k) square covariance ever formed.
        const arma::vec rz = r * z;
        const double q = arma::dot(z, rz) + 1.0;
        const arma::vec prediction = m.t() * z;
        const arma::vec e = y - prediction;

        if (mean != nullptr)
        {
            *mean = prediction;
        }
        const double score = log_mvt(e, q * (d_pred / df_pred), df_pred);

        const arma::vec a = rz / q;
        m += a * e.t();
        c = arma::symmatu(r - (a * a.t()) * q);
        df = df_pred + 1.0;
        d = arma::symmatu(d_pred + (e * e.t()) / q);

        return score;
    }
};

/// What the filter and its retrospective pass leave behind, in the flattened
/// per-period layout both posterior structs store: one column per period
/// throughout, as a draw is one column elsewhere.
struct DiscountPath
{
    arma::mat a;       ///< (n_design k) x tt, the posterior mean of vec(B_t).
    arma::mat a_scale; ///< (n_design k) x tt, the marginal Student t scale.
    arma::mat a_cov;   ///< (n_design n_design) x tt, the regressor side of it.
    arma::mat u_sigma; ///< (k k) x tt, the posterior estimate of Sigma_t.
    arma::vec df;      ///< tt.

    arma::mat forecast_mean; ///< k x tt, the one step ahead predictive mean.
    arma::vec loglik;        ///< tt, its exact log density at what happened.
};

/// Runs the filter over `y`, tt x k, against the design `x`, tt x n_design, and
/// smooths it.
///
/// Reports progress once per period and honours an interrupt thrown from the
/// reporter. Consumes no random numbers at all, which is what lets a host run
/// it twice and get the same bits.
inline DiscountPath run_discount_filter(const arma::mat &y, const arma::mat &x,
                                        const MatrixNormalPrior &a_prior,
                                        const WishartPrior &u_sigma_prior, const double d_beta,
                                        const double d_sigma, Reporter &reporter)
{
    const arma::uword tt = y.n_rows;
    const arma::uword k = y.n_cols;
    const arma::uword n_design = x.n_cols;

    DiscountState state;
    state.m = a_prior.mean;
    state.c = a_prior.cov;
    state.df = static_cast<double>(u_sigma_prior.df);
    state.d = state.df * u_sigma_prior.scale;

    arma::cube m_path(n_design, k, tt, arma::fill::zeros);
    arma::cube c_path(n_design, n_design, tt, arma::fill::zeros);
    arma::cube s_path(k, k, tt, arma::fill::zeros);
    arma::vec df_path(tt, arma::fill::zeros);

    DiscountPath out;
    out.forecast_mean.set_size(k, tt);
    out.loglik.set_size(tt);

    arma::vec mean(k);
    for (arma::uword t = 0; t < tt; t++)
    {
        reporter.check_interrupt();

        out.loglik(t) = state.step(x.row(t).t(), y.row(t).t(), d_beta, d_sigma, &mean);
        out.forecast_mean.col(t) = mean;

        m_path.slice(t) = state.m;
        c_path.slice(t) = state.c;
        s_path.slice(t) = state.d / state.df;
        df_path(t) = state.df;

        reporter.progress(static_cast<long long>(t) + 1, static_cast<long long>(tt));
    }

    // The retrospective pass. For a random walk whose predicted covariance is
    // C_{t-1} / delta the Rauch-Tung-Striebel gain is exactly delta times the
    // identity, so smoothing the mean is backward exponential smoothing and the
    // covariance recursion collapses with it. The precision smooths
    // harmonically, which is the retrospective distribution of a discounted
    // Wishart.
    arma::cube m_s = m_path;
    arma::cube c_s = c_path;
    arma::cube s_s = s_path;
    arma::vec df_s = df_path;

    for (arma::uword t = tt - 1; t-- > 0;)
    {
        reporter.check_interrupt();

        m_s.slice(t) = (1.0 - d_beta) * m_path.slice(t) + d_beta * m_s.slice(t + 1);
        c_s.slice(t) = arma::symmatu((1.0 - d_beta) * c_path.slice(t) +
                                     d_beta * d_beta * c_s.slice(t + 1));

        df_s(t) = (1.0 - d_sigma) * df_path(t) + d_sigma * df_s(t + 1);
        const arma::mat precision = (1.0 - d_sigma) * arma::inv_sympd(s_path.slice(t)) +
                                    d_sigma * arma::inv_sympd(s_s.slice(t + 1));
        s_s.slice(t) = arma::symmatu(arma::inv_sympd(arma::symmatu(precision)));
    }

    out.a.set_size(n_design * k, tt);
    out.a_scale.set_size(n_design * k, tt);
    out.a_cov.set_size(n_design * n_design, tt);
    out.u_sigma.set_size(k * k, tt);
    out.df = df_s;

    for (arma::uword t = 0; t < tt; t++)
    {
        out.a.col(t) = arma::vectorise(m_s.slice(t).t());
        out.a_scale.col(t) = coefficient_scale(c_s.slice(t), s_s.slice(t));
        out.a_cov.col(t) = arma::vectorise(c_s.slice(t));
        out.u_sigma.col(t) = arma::vectorise(s_s.slice(t));
    }

    return out;
}

/// I.i.d. draws from the posterior of one period: first `Sigma ~ IW(df, df S)`,
/// then `vec(Theta) | Sigma ~ N(m, Sigma kron C)`. Returns (n_design k) x draws,
/// the layout every sampler in this project returns.
///
/// Correct for that period alone. The smoothed posterior is not independent
/// across periods, so calling this per period and joining the columns does not
/// give a draw of the coefficient path.
inline arma::mat draw_discount_period(const arma::mat &a, const arma::mat &a_cov,
                                      const arma::mat &u_sigma, const arma::vec &df_path,
                                      const arma::uword period, const arma::uword draws,
                                      const arma::uword k)
{
    if (period >= u_sigma.n_cols)
    {
        throw std::invalid_argument("period is past the end of the posterior");
    }
    if (draws == 0)
    {
        throw std::invalid_argument("a draw count of zero draws nothing");
    }

    const arma::uword n_design = a.n_rows / k;
    const double df = df_path(period);

    const arma::mat s = arma::reshape(u_sigma.col(period), k, k);
    const arma::mat c = arma::reshape(a_cov.col(period), n_design, n_design);
    const arma::mat m = arma::reshape(a.col(period), k, n_design).t();

    arma::mat chol_inv_scale;
    if (!arma::chol(chol_inv_scale, arma::inv_sympd(arma::symmatu(df * s)), "lower"))
    {
        throw std::invalid_argument("the posterior scale of this period is not positive definite");
    }
    arma::mat chol_c;
    if (!arma::chol(chol_c, arma::symmatu(c), "lower"))
    {
        throw std::invalid_argument(
            "the coefficient covariance of this period is not positive definite");
    }

    arma::mat out(n_design * k, draws);
    for (arma::uword i = 0; i < draws; i++)
    {
        const arma::mat sigma = draw_inverse_wishart(chol_inv_scale, df, k);
        arma::mat chol_sigma;
        if (!arma::chol(chol_sigma, arma::symmatu(sigma), "lower"))
        {
            throw std::invalid_argument("a drawn error covariance was not positive definite");
        }
        // vec(Theta) ~ N(vec(m), Sigma kron C): one standard normal matrix
        // scaled on the left by C's factor and on the right by Sigma's.
        const arma::mat noise(n_design, k, arma::fill::randn);
        const arma::mat theta = m + chol_c * noise * chol_sigma.t();
        out.col(i) = arma::vectorise(theta.t());
    }
    return out;
}

/// The lower Cholesky factor of the random walk innovation the coefficient
/// discount implies at a period whose posterior covariance is `c`, and zero at
/// a discount of one, where the walk does not move.
///
/// A factorisation that fails is not an error here: it means a direction the
/// sample pinned down to working precision, and a forecast that does not move
/// the coefficients along it is the right answer rather than a fallback.
inline arma::mat discount_walk_root(const arma::mat &c, const double d_beta)
{
    arma::mat root;
    if (!arma::chol(root, arma::symmatu(((1.0 - d_beta) / d_beta) * c), "lower"))
    {
        root.zeros(c.n_rows, c.n_cols);
    }
    return root;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_DISCOUNT_SUPPORT_H
