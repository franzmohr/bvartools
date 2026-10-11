// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_ADAPTIVE_WALK_H
#define BAYESTS_CORE_MODELS_ADAPTIVE_WALK_H

#include "bayests/arma.h"

#include <algorithm>
#include <cmath>

namespace bayests::core
{

/// The proposal of a random walk Metropolis-Hastings step, adapted over the
/// burn-in and held after it, so that the kept draws come from one fixed
/// kernel. Serves UcNormalSkew's scale coefficients and the decay of
/// VarNormalWishart's volatility break.
///
/// The adaptation is Andrieu and Thoms (2008, algorithm 4): the proposal is
/// N(theta, s^2 C), with C the running covariance of the chain -- Haario,
/// Saksman and Tamminen (2001) -- once there are enough draws to estimate it and
/// the identity before, and log s moved after every proposal by a Robbins-Monro
/// step towards an acceptance rate of 0.234. The scale is what keeps the first
/// half of that working: a covariance learnt from a chain that never moved is
/// the zero matrix, and the identity is too wide for a long sample, so a
/// covariance alone stalls exactly where it is needed.
class AdaptiveWalk
{
public:
    explicit AdaptiveWalk(const arma::uword d)
        : d_(d), mean_(d, arma::fill::zeros), scatter_(d, d, arma::fill::zeros),
          shape_(arma::eye<arma::mat>(d, d)), log_scale_(std::log(0.1))
    {
    }

    arma::vec propose(const arma::vec &theta) const
    {
        return theta + std::exp(log_scale_) * (shape_ * arma::randn<arma::vec>(d_));
    }

    /// One burn-in step: the scale towards the target acceptance, then the
    /// chain's moments, and from 10 draws a dimension on the proposal's shape.
    void learn(const arma::vec &theta, const double acceptance)
    {
        n_++;
        const double gain = std::pow(static_cast<double>(n_), -0.6);
        log_scale_ += gain * (std::min(acceptance, 1.0) - 0.234);

        const arma::vec delta = theta - mean_;
        mean_ += delta / static_cast<double>(n_);
        scatter_ += delta * arma::trans(theta - mean_);
        if (n_ >= 10 * d_)
        {
            const arma::mat cov = scatter_ / static_cast<double>(n_ - 1);
            arma::mat root;
            if (arma::chol(root, arma::symmatu(cov + arma::eye<arma::mat>(d_, d_) * 1e-10),
                           "lower"))
            {
                // The scale was tuned for the identity; for the chain's own
                // covariance the optimal one is near 2.38 / sqrt(d) (Roberts,
                // Gelman and Gilks 1997), so it starts again from there.
                if (!learnt_)
                {
                    log_scale_ = std::log(2.38 / std::sqrt(static_cast<double>(d_)));
                    learnt_ = true;
                }
                shape_ = root;
            }
        }
    }

private:
    arma::uword d_ = 0;
    arma::uword n_ = 0;
    arma::vec mean_;
    arma::mat scatter_;
    arma::mat shape_;
    double log_scale_ = 0.0;
    bool learnt_ = false;
};

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_ADAPTIVE_WALK_H
