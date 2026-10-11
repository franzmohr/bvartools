// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_VOLATILITY_BREAK_SUPPORT_H
#define BAYESTS_CORE_MODELS_VOLATILITY_BREAK_SUPPORT_H

#include "bayests/arma.h"
#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/reporter.h"
#include "bayests/spec.h"
#include "core/models/adaptive_walk.h"
#include "core/models/error_scale_checks.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace bayests::core
{

// The volatility break of Lenza and Primiceri (2022): from a known period the
// error covariance is s_t^2 Sigma, a common scale on every shock. See
// VarSpec::volatility_break for the path and VolatilityBreakPrior for the prior.
//
// Conditional on the scales the model is the VAR of the rescaled data,
// y_t / s_t on x_t / s_t, so the coefficient and Wishart blocks are the ones a
// sampler already has, run on rescaled rows. What this file adds is the path
// itself, the draw of its parameters given the residuals, and the refusals.

/// s_t for the periods first .. first + n - 1, counted from zero along
/// /data/train/y and on past its end into a forecast horizon. `theta` holds the
/// volatility_scales free scales, then the decay rate. All ones without a break.
inline arma::vec volatility_path(const VarSpec &spec, const arma::vec &theta, const arma::uword first,
                                 const arma::uword n)
{
    arma::vec s(n, arma::fill::ones);
    if (!spec.uses_volatility_break())
    {
        return s;
    }
    const long long start = spec.volatility_break - 1;
    const long long m = spec.volatility_scales;
    const double last = theta(static_cast<arma::uword>(m - 1));
    const double rho = theta(static_cast<arma::uword>(m));
    for (arma::uword i = 0; i < n; i++)
    {
        const long long j = static_cast<long long>(first + i) - start;
        if (j < 0)
        {
            continue;
        }
        s(i) = j < m ? theta(static_cast<arma::uword>(j))
                     : 1.0 + (last - 1.0) * std::pow(rho, static_cast<double>(j - m + 1));
    }
    return s;
}

/// One draw from Gamma(shape a, rate b) restricted to (0, 1].
///
/// Two exact rejection samplers, and the one that accepts more often: the
/// unrestricted gamma kept when it lands in (0, 1], accepting with
/// probability P(a, b) = gamma(a, b) / Gamma(a); or w = U^(1/a), whose density
/// a w^(a-1) the target's is e^(-b w) times, accepting with probability
/// a gamma(a, b) / b^a. The ratio of the second to the first is
/// Gamma(a + 1) / b^a whatever the incomplete gamma is, so the choice needs
/// no special function.
inline double draw_gamma_below_one(const double a, const double b)
{
    const bool from_power = a * std::log(b) < std::lgamma(a + 1.0);
    for (long tries = 0; tries < 100000000L; tries++)
    {
        if (from_power)
        {
            const double w = std::pow(arma::randu<double>(), 1.0 / a);
            if (arma::randu<double>() <= std::exp(-b * w))
            {
                return w;
            }
        }
        else
        {
            const double w = arma::randg<double>(arma::distr_param(a, 1.0 / b));
            if (w <= 1.0)
            {
                return w;
            }
        }
    }
    throw std::runtime_error("no draw of a gamma restricted to (0, 1] was accepted in 1e8 tries, "
                             "shape " + std::to_string(a) + " and rate " + std::to_string(b));
}

/// The volatility block of a sampler: the path s_t over the sample and the
/// draw of its parameters given the residuals and the error precision.
///
/// Given the rest, each free scale but the last moves one period only, and its
/// conditional is exact: with q_t = u_t' Sigma^-1 u_t, s_t^2 is inverse gamma
/// with shape (k + shape) / 2 and scale q_t / 2, restricted to at least one,
/// drawn as the reciprocal of a gamma restricted to (0, 1]. A free scale whose
/// period lies past the sample is drawn from its Pareto prior. The last free
/// scale and the decay rate move every period from there on together and have
/// no closed form; they are one random walk Metropolis-Hastings step in
/// (log(s - 1), logit rho), its proposal adapted over the burn-in
/// (AdaptiveWalk) and held after it.
///
/// Inactive without a break: it holds ones and draws nothing, so a sampler
/// that owns one consumes no random number for it.
class VolatilityBreak
{
public:
    VolatilityBreak(const VarSpec &spec, const VolatilityBreakPrior &prior, const arma::vec &initial,
                    const arma::uword tt, const arma::uword k)
        : spec_(spec), prior_(prior), tt_(tt), k_(k), walk_(2), scales_(tt, arma::fill::ones)
    {
        if (!spec.uses_volatility_break())
        {
            return;
        }
        active_ = true;
        m_ = static_cast<arma::uword>(spec.volatility_scales);
        start_ = static_cast<arma::uword>(spec.volatility_break - 1);
        if (initial.n_elem > 0)
        {
            theta_ = initial;
        }
        else
        {
            theta_ = arma::vec(m_ + 1);
            theta_.head(m_).fill(2.0);
            theta_(m_) = 0.5;
        }
        scales_ = volatility_path(spec_, theta_, 0, tt_);
    }

    bool active() const { return active_; }

    /// s_t over the sample, tt of them.
    const arma::vec &scales() const { return scales_; }

    /// The free scales, then the decay rate.
    const arma::vec &theta() const { return theta_; }

    /// The first period that is scaled, counted from zero.
    arma::uword start() const { return start_; }

    /// One draw of the parameters given `u`, the k x tt residuals on the
    /// original scale, and the error precision. `burnin` says whether the
    /// proposal is still being adapted.
    void update(const arma::mat &u, const arma::mat &u_sigma_inv, const bool burnin)
    {
        if (!active_)
        {
            return;
        }
        arma::vec q(tt_, arma::fill::zeros);
        const arma::mat tail = u.cols(start_, tt_ - 1);
        q.subvec(start_, tt_ - 1) = arma::trans(arma::sum((u_sigma_inv * tail) % tail, 0));

        const double kd = static_cast<double>(k_);
        for (arma::uword j = 0; j + 1 < m_; j++)
        {
            const arma::uword t = start_ + j;
            theta_(j) = t < tt_
                            ? 1.0 / std::sqrt(draw_gamma_below_one((kd + prior_.shape) / 2.0, q(t) / 2.0))
                            : std::pow(arma::randu<double>(), -1.0 / prior_.shape);
        }

        const arma::vec current = {std::log(theta_(m_ - 1) - 1.0),
                                   std::log(theta_(m_) / (1.0 - theta_(m_)))};
        const arma::vec proposal = walk_.propose(current);
        const double log_ratio = log_target(proposal, q) - log_target(current, q);
        const bool accept = std::log(arma::randu<double>()) < log_ratio;
        const arma::vec &kept = accept ? proposal : current;
        theta_(m_ - 1) = 1.0 + std::exp(kept(0));
        theta_(m_) = 1.0 / (1.0 + std::exp(-kept(1)));
        if (burnin)
        {
            walk_.learn(kept, std::exp(std::min(log_ratio, 0.0)));
        }
        else
        {
            proposals_++;
            accepted_ += accept ? 1 : 0;
        }
        scales_ = volatility_path(spec_, theta_, 0, tt_);
    }

    /// What the step did, through the reporter: its acceptance after the
    /// burn-in, and a burn-in too short for the proposal to adapt.
    void report(Reporter &reporter, const int burnin) const
    {
        if (!active_)
        {
            return;
        }
        if (proposals_ > 0)
        {
            reporter.message("the volatility break's Metropolis-Hastings step accepted " +
                             std::to_string(accepted_) + " of " + std::to_string(proposals_) +
                             " proposals after the burn-in");
        }
        if (burnin < 20)
        {
            reporter.message("a burn-in of " + std::to_string(burnin) +
                             " is too short for the volatility break's proposal to adapt, which "
                             "takes 20 draws; it stayed a scaled identity");
        }
    }

private:
    /// The log conditional of (log(s_{m-1} - 1), logit rho), Jacobian included:
    /// every period from the last free scale's on, then the Pareto and the
    /// beta.
    double log_target(const arma::vec &phi, const arma::vec &q) const
    {
        const double last = 1.0 + std::exp(phi(0));
        const double rho = 1.0 / (1.0 + std::exp(-phi(1)));
        const double kd = static_cast<double>(k_);
        double out = 0.0;
        for (arma::uword t = start_ + m_ - 1; t < tt_; t++)
        {
            const double j = static_cast<double>(t - start_) - static_cast<double>(m_ - 1);
            const double s = 1.0 + (last - 1.0) * std::pow(rho, j);
            out += -kd * std::log(s) - q(t) / (2.0 * s * s);
        }
        out += -(prior_.shape + 1.0) * std::log(last) + phi(0);
        out += prior_.rho_shape1 * std::log(rho) + prior_.rho_shape2 * std::log(1.0 - rho);
        return out;
    }

    const VarSpec &spec_;
    const VolatilityBreakPrior &prior_;
    arma::uword tt_ = 0;
    arma::uword k_ = 0;
    bool active_ = false;
    arma::uword m_ = 0;
    arma::uword start_ = 0;
    AdaptiveWalk walk_;
    arma::vec theta_;
    arma::vec scales_;
    long long proposals_ = 0;
    long long accepted_ = 0;
};

/// Scales the rows of a stacked response and its SUR regressors -- k rows per
/// period, period-major -- by 1 / s_t, from the break on. `y` and `z` are the
/// original rows; the outputs are overwritten.
inline void rescale_rows(const arma::vec &y, const arma::mat &z, const arma::vec &s, const arma::uword k,
                         const arma::uword start, arma::vec &y_out, arma::mat &z_out)
{
    const arma::uword tt = s.n_elem;
    for (arma::uword t = start; t < tt; t++)
    {
        const arma::uword r0 = t * k;
        const arma::uword r1 = r0 + k - 1;
        y_out.subvec(r0, r1) = y.subvec(r0, r1) / s(t);
        if (z.n_cols > 0)
        {
            z_out.rows(r0, r1) = z.rows(r0, r1) / s(t);
        }
    }
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_VOLATILITY_BREAK_SUPPORT_H
