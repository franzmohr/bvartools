// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_STUDENT_T_SUPPORT_H
#define BAYESTS_CORE_MODELS_STUDENT_T_SUPPORT_H

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

// Student t errors as a scale mixture: u_t ~ N(0, lambda_t Sigma) with
// lambda_t ~ IG(nu / 2, nu / 2) independently, so that u_t is multivariate t
// with nu degrees of freedom and scale Sigma. The scale is common to every
// equation, which is what keeps the Wishart block conjugate -- given lambda the
// model is the VAR of the rows divided by sqrt(lambda_t), the same rescaling a
// volatility break uses, and the two multiply: lambda_t s_t^2 Sigma. A period
// the shocks of which are large gets a large lambda_t and so a small weight,
// without anyone having to say which periods those are.

/// The log density of the k-dimensional Student t with `nu` degrees of
/// freedom and scale matrix S at a point whose Mahalanobis distance under S is
/// `q`; `log_det_s` is log |S|.
inline double student_t_log_density(const double q, const double log_det_s, const double k,
                                    const double nu)
{
    return std::lgamma(0.5 * (nu + k)) - std::lgamma(0.5 * nu) -
           0.5 * k * std::log(nu * arma::datum::pi) - 0.5 * log_det_s -
           0.5 * (nu + k) * std::log1p(q / nu);
}

/// The Student t block of a sampler: the mixing scales lambda_t over the
/// sample and the degrees of freedom.
///
/// Given everything else, lambda_t is inverse gamma with shape (nu + k) / 2 and
/// scale (nu + q_t) / 2, q_t = e_t' Sigma^-1 e_t for the residual e_t of period
/// t with any volatility break already divided out -- an exact draw. Given the
/// lambdas, nu has no closed form; under a gamma prior it is one random walk
/// Metropolis-Hastings step in log nu, adapted over the burn-in (AdaptiveWalk)
/// and held after it. Under a fixed nu the step is skipped.
///
/// Inactive without /model/student_t: it holds ones and draws nothing, so a
/// sampler that owns one consumes no random number for it.
class StudentTErrors
{
public:
    StudentTErrors(const VarSpec &spec, const StudentTPrior &prior, const double initial_df,
                   const arma::uword tt, const arma::uword k)
        : prior_(prior), tt_(tt), k_(k), walk_(1), lambda_(tt, arma::fill::ones)
    {
        if (!spec.student_t)
        {
            return;
        }
        active_ = true;
        fixed_ = prior.df != 0.0;
        if (fixed_)
        {
            df_ = prior.df;
        }
        else
        {
            df_ = initial_df != 0.0 ? initial_df : prior.shape / prior.rate;
        }
    }

    bool active() const { return active_; }

    /// lambda_t over the sample, tt of them.
    const arma::vec &lambda() const { return lambda_; }

    /// The degrees of freedom.
    double df() const { return df_; }

    /// One draw given `e`, the k x tt residuals with any volatility break
    /// divided out, and the error precision.
    void update(const arma::mat &e, const arma::mat &u_sigma_inv, const bool burnin)
    {
        if (!active_)
        {
            return;
        }
        const arma::rowvec q = arma::sum((u_sigma_inv * e) % e, 0);
        const double shape = 0.5 * (df_ + static_cast<double>(k_));
        for (arma::uword t = 0; t < tt_; t++)
        {
            const double rate = 0.5 * (df_ + q(t));
            lambda_(t) = 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / rate));
        }
        if (fixed_)
        {
            return;
        }

        const arma::vec current = {std::log(df_)};
        const arma::vec proposal = walk_.propose(current);
        const double log_ratio = log_target(proposal(0)) - log_target(current(0));
        const bool accept = std::log(arma::randu<double>()) < log_ratio;
        const arma::vec &kept = accept ? proposal : current;
        df_ = std::exp(kept(0));
        if (burnin)
        {
            walk_.learn(kept, std::exp(std::min(log_ratio, 0.0)));
        }
        else
        {
            proposals_++;
            accepted_ += accept ? 1 : 0;
        }
    }

    /// What the step did, through the reporter.
    void report(Reporter &reporter, const int burnin) const
    {
        if (!active_ || fixed_)
        {
            return;
        }
        if (proposals_ > 0)
        {
            reporter.message("the Student t degrees of freedom's Metropolis-Hastings step "
                             "accepted " +
                             std::to_string(accepted_) + " of " + std::to_string(proposals_) +
                             " proposals after the burn-in");
        }
        if (burnin < 10)
        {
            reporter.message("a burn-in of " + std::to_string(burnin) +
                             " is too short for the Student t degrees of freedom's proposal to "
                             "adapt, which takes 10 draws; it stayed a scaled identity");
        }
    }

private:
    /// log p(nu | lambda) in log nu, Jacobian included: every lambda_t's
    /// inverse gamma with shape and scale nu / 2, then the gamma prior.
    double log_target(const double log_nu) const
    {
        const double nu = std::exp(log_nu);
        const double half = 0.5 * nu;
        double out = static_cast<double>(tt_) * (half * std::log(half) - std::lgamma(half));
        out -= arma::accu((half + 1.0) * arma::log(lambda_) + half / lambda_);
        out += (prior_.shape - 1.0) * std::log(nu) - prior_.rate * nu + log_nu;
        return out;
    }

    const StudentTPrior &prior_;
    arma::uword tt_ = 0;
    arma::uword k_ = 0;
    bool active_ = false;
    bool fixed_ = false;
    double df_ = 0.0;
    AdaptiveWalk walk_;
    arma::vec lambda_;
    long long proposals_ = 0;
    long long accepted_ = 0;
};

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_STUDENT_T_SUPPORT_H
