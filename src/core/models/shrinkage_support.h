// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_SHRINKAGE_SUPPORT_H
#define BAYESTS_CORE_MODELS_SHRINKAGE_SUPPORT_H

#include "bayests/priors.h"
#include "bayests/reporter.h"
#include "bayests/spec.h"
#include "core/algorithms/gig_hormann_leydold_2014.h"
#include "core/models/model_support.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace bayests::core
{

// Adaptive priors on a block of coefficients, and the stationarity condition,
// for the constant-coefficient Gaussian VARs. Both act on the coefficient draw
// alone: shrinkage changes the prior it is made under, stationarity which draws
// are kept. Neither consumes a random number where it is not asked for, which
// is what leaves every other file's draws where they were.

/// An inverse gamma draw, IG(shape, rate): the reciprocal of a gamma draw of
/// that shape and rate.
inline double draw_inverse_gamma(const double shape, const double rate)
{
    return 1.0 / arma::randg<double>(arma::distr_param(shape, 1.0 / rate));
}

/// Whether the scheme keeps a local scale per coefficient beside its group
/// scales, which the samplers then write as /posterior/<block>/local.
inline bool shrinkage_has_local(const Shrinkage type)
{
    return type == Shrinkage::horseshoe || type == Shrinkage::normal_gamma;
}

/// Floors of the normal-gamma draw. The squared deviation a local variance is
/// drawn from is floored so that its generalised inverse Gaussian stays proper
/// when theta < 1/2 and a coefficient lands on its prior mean; the local
/// variance itself so that the prior precision, the base precision divided by
/// it, stays finite. A coefficient whose prior standard deviation is 1e-6 of
/// its base one is zero for every purpose a VAR has, so the second floor moves
/// nothing but a division by zero.
constexpr double kNormalGammaMinDeviation = 1e-300;
constexpr double kNormalGammaMinLocal = 1e-12;

/// The Metropolis-Hastings step on log theta of a normal-gamma prior that
/// draws theta: its starting step size, and the acceptance rate the burn-in
/// tunes it towards.
constexpr double kThetaInitialStep = 0.5;
constexpr double kThetaTargetAcceptance = 0.3;

/// The adaptive prior of VarSpec::shrinkage over one block of coefficients,
/// and the scales it keeps from one sweep to the next.
///
/// Before every draw of the block the sampler asks for precision() and rhs(),
/// the prior it is to be drawn under; after it, update() draws the scales given
/// the block. Inactive (Shrinkage::none) it holds nothing and the sampler keeps
/// the prior it had.
///
/// Under a normal-gamma prior that draws theta, the first `adapt_sweeps` calls
/// of update() also tune the step of the Metropolis-Hastings move on log theta
/// towards an acceptance rate of kThetaTargetAcceptance; the sampler passes its
/// burn-in, so the step is fixed over every kept draw and the chain there is an
/// ordinary Metropolis-within-Gibbs one.
class CoefficientShrinkage
{
public:
    CoefficientShrinkage(const Shrinkage type, const ShrinkagePrior &prior, const NormalPrior &base,
                         const arma::vec &initial_scale, const arma::vec &initial_local,
                         const arma::vec &initial_theta = arma::vec(), const int adapt_sweeps = 0)
        : type_(type), prior_(prior), base_(base), adapt_sweeps_(adapt_sweeps)
    {
        if (type_ == Shrinkage::none)
        {
            return;
        }
        groups_ = prior_.group.max();
        scale_ = initial_scale.n_elem > 0 ? initial_scale : arma::vec(groups_, arma::fill::ones);
        if (shrinkage_has_local(type_))
        {
            local_ = initial_local.n_elem > 0 ? initial_local
                                              : arma::vec(prior_.group.n_elem, arma::fill::ones);
        }
        if (type_ == Shrinkage::horseshoe)
        {
            local_aux_.ones(prior_.group.n_elem);
            global_aux_.ones(groups_);
        }
        if (type_ == Shrinkage::normal_gamma)
        {
            if (draws_theta())
            {
                theta_ = initial_theta.n_elem > 0 ? initial_theta : arma::vec(1.0 / prior_.theta_rate);
                log_step_ = arma::vec(groups_, arma::fill::value(std::log(kThetaInitialStep)));
                accepted_.zeros(groups_);
            }
            else
            {
                theta_ = prior_.theta;
            }
        }
        base_diag_ = base_.v_inv.diag();
    }

    bool active() const { return type_ != Shrinkage::none; }
    arma::uword groups() const { return groups_; }

    /// The group scales: s_g under `minnesota`, tau_g^2 under `horseshoe`,
    /// the global rate lambda_g under `normal_gamma`.
    const arma::vec &scale() const { return scale_; }

    /// The local scales lambda_j^2 under `horseshoe` and the local variances
    /// psi_j under `normal_gamma`, one per coefficient, and one for an unshrunk
    /// coefficient too, left at its start.
    const arma::vec &local() const { return local_; }

    /// Whether theta_g of a normal-gamma prior is drawn rather than fixed.
    bool draws_theta() const
    {
        return type_ == Shrinkage::normal_gamma && prior_.theta_rate.n_elem > 0;
    }

    /// theta_g under `normal_gamma`, fixed or the current draw.
    const arma::vec &theta() const { return theta_; }

    /// The step of the move on log theta_g, and the share of moves accepted
    /// since the burn-in ended; for tests and diagnostics.
    arma::vec theta_step() const { return arma::exp(log_step_); }
    arma::vec theta_acceptance() const
    {
        const int kept = sweeps_ - adapt_sweeps_;
        return kept > 0 ? arma::vec(accepted_ / kept) : arma::vec(groups_, arma::fill::zeros);
    }

    /// The prior precision of the block the next draw is made under: the base
    /// precision, its diagonal divided by the current prior variance multiplier
    /// wherever a coefficient is shrunk.
    arma::mat precision() const
    {
        arma::mat out = base_.v_inv;
        for (arma::uword j = 0; j < prior_.group.n_elem; j++)
        {
            if (prior_.group(j) > 0)
            {
                out(j, j) = base_diag_(j) / multiplier(j);
            }
        }
        return out;
    }

    /// precision() times the prior mean: the prior's part of the right hand side.
    arma::vec rhs() const { return precision() * base_.mu; }

    /// One Gibbs draw of every scale given the block `a`.
    void update(const arma::vec &a)
    {
        if (type_ == Shrinkage::none)
        {
            return;
        }
        // (a_j - mu_j)^2 v_inv(j, j): the squared deviation in units of the base
        // prior variance, what every conditional below is made of.
        arma::vec dev(prior_.group.n_elem);
        for (arma::uword j = 0; j < dev.n_elem; j++)
        {
            const double d = a(j) - base_.mu(j);
            dev(j) = d * d * base_diag_(j);
        }

        if (type_ == Shrinkage::normal_gamma)
        {
            // psi_j | a_j, lambda_g ~ GIG(theta_g - 1/2, dev_j, theta_g lambda_g),
            // then lambda_g | psi ~ G(shape_g + theta_g n_g,
            // rate_g + theta_g sum_j psi_j / 2), after Huber and Feldkircher
            // (2019).
            arma::vec count(groups_, arma::fill::zeros), sum(groups_, arma::fill::zeros);
            for (arma::uword j = 0; j < dev.n_elem; j++)
            {
                const arma::uword g = prior_.group(j);
                if (g == 0)
                {
                    continue;
                }
                const double theta = theta_(g - 1);
                local_(j) = std::max(
                    gig_hormann_leydold_2014(theta - 0.5, std::max(dev(j), kNormalGammaMinDeviation),
                                             theta * scale_(g - 1)),
                    kNormalGammaMinLocal);
                count(g - 1) += 1.0;
                sum(g - 1) += local_(j);
            }
            for (arma::uword g = 0; g < groups_; g++)
            {
                const double theta = theta_(g);
                scale_(g) = arma::randg<double>(arma::distr_param(
                    prior_.shape(g) + theta * count(g), 1.0 / (prior_.rate(g) + 0.5 * theta * sum(g))));
            }
            if (draws_theta())
            {
                update_theta();
            }
            sweeps_++;
            return;
        }

        if (type_ == Shrinkage::minnesota)
        {
            arma::vec count(groups_, arma::fill::zeros), sum(groups_, arma::fill::zeros);
            for (arma::uword j = 0; j < dev.n_elem; j++)
            {
                if (prior_.group(j) > 0)
                {
                    count(prior_.group(j) - 1) += 1.0;
                    sum(prior_.group(j) - 1) += dev(j);
                }
            }
            for (arma::uword g = 0; g < groups_; g++)
            {
                scale_(g) = draw_inverse_gamma(prior_.shape(g) + 0.5 * count(g),
                                               prior_.rate(g) + 0.5 * sum(g));
            }
            return;
        }

        // Horseshoe, after Makalic and Schmidt (2016): with nu_j and xi_g the
        // auxiliary variables of the half-Cauchy scales, every conditional is
        // an inverse gamma.
        arma::vec count(groups_, arma::fill::zeros), sum(groups_, arma::fill::zeros);
        for (arma::uword j = 0; j < dev.n_elem; j++)
        {
            const arma::uword g = prior_.group(j);
            if (g == 0)
            {
                continue;
            }
            local_(j) = draw_inverse_gamma(1.0, 1.0 / local_aux_(j) + 0.5 * dev(j) / scale_(g - 1));
            local_aux_(j) = draw_inverse_gamma(1.0, 1.0 + 1.0 / local_(j));
            count(g - 1) += 1.0;
            sum(g - 1) += dev(j) / local_(j);
        }
        for (arma::uword g = 0; g < groups_; g++)
        {
            scale_(g) = draw_inverse_gamma(0.5 * (count(g) + 1.0), 1.0 / global_aux_(g) + 0.5 * sum(g));
            global_aux_(g) = draw_inverse_gamma(1.0, 1.0 + 1.0 / scale_(g));
        }
    }

private:
    /// log p(theta_g | psi, lambda_g) up to a constant: the exponential prior
    /// times the gamma density of every local variance of the group,
    ///
    ///     -r theta + n (theta log(theta lambda / 2) - lgamma(theta))
    ///       + (theta - 1) sum log psi_j - theta lambda sum psi_j / 2.
    double theta_log_density(const double theta, const arma::uword g, const double n,
                             const double sum_log, const double sum) const
    {
        const double lambda = scale_(g);
        return -prior_.theta_rate(g) * theta +
               n * (theta * std::log(0.5 * theta * lambda) - std::lgamma(theta)) +
               (theta - 1.0) * sum_log - 0.5 * theta * lambda * sum;
    }

    /// One random walk Metropolis-Hastings move on log theta_g per group; the
    /// Jacobian of the log is the log theta term of the acceptance ratio.
    /// During the first adapt_sweeps_ sweeps the step follows a Robbins-Monro
    /// recursion towards kThetaTargetAcceptance; after them it is fixed and the
    /// acceptances are counted.
    void update_theta()
    {
        arma::vec n(groups_, arma::fill::zeros), sum_log(groups_, arma::fill::zeros),
            sum(groups_, arma::fill::zeros);
        for (arma::uword j = 0; j < prior_.group.n_elem; j++)
        {
            const arma::uword g = prior_.group(j);
            if (g > 0)
            {
                n(g - 1) += 1.0;
                sum_log(g - 1) += std::log(local_(j));
                sum(g - 1) += local_(j);
            }
        }
        const bool adapting = sweeps_ < adapt_sweeps_;
        for (arma::uword g = 0; g < groups_; g++)
        {
            const double current = theta_(g);
            const double proposal = current * std::exp(std::exp(log_step_(g)) * arma::randn<double>());
            const double log_ratio =
                theta_log_density(proposal, g, n(g), sum_log(g), sum(g)) -
                theta_log_density(current, g, n(g), sum_log(g), sum(g)) + std::log(proposal / current);
            const double u = arma::randu<double>();
            const bool accept = std::log(u) < log_ratio;
            if (accept)
            {
                theta_(g) = proposal;
            }
            if (adapting)
            {
                log_step_(g) += ((accept ? 1.0 : 0.0) - kThetaTargetAcceptance) /
                                std::sqrt(1.0 + static_cast<double>(sweeps_));
            }
            else if (accept)
            {
                accepted_(g) += 1.0;
            }
        }
    }

    double multiplier(const arma::uword j) const
    {
        const arma::uword g = prior_.group(j) - 1;
        switch (type_)
        {
        case Shrinkage::minnesota:
            return scale_(g);
        case Shrinkage::normal_gamma:
            return local_(j);
        default:
            return scale_(g) * local_(j);
        }
    }

    Shrinkage type_;
    const ShrinkagePrior &prior_;
    const NormalPrior &base_;
    arma::uword groups_ = 0;
    arma::vec base_diag_;
    arma::vec scale_;
    arma::vec local_;
    arma::vec local_aux_;
    arma::vec global_aux_;
    arma::vec theta_;
    arma::vec log_step_;
    arma::vec accepted_;
    int adapt_sweeps_ = 0;
    int sweeps_ = 0;
};

/// Whether VarSpec::shrinkage and VarSpec::stationary are read by the
/// algorithm: refused where not, rather than ignored, and refused together
/// with what rearranges the prior or the draw they act on.
inline void require_supported_shrinkage(const VarSpec &spec, const bool supported,
                                        const char *algorithm)
{
    const bool asked = spec.shrinkage != Shrinkage::none || spec.stationary;
    if (!asked)
    {
        return;
    }
    const std::string name(algorithm);
    if (!supported)
    {
        throw std::invalid_argument(
            name + " does not read /model/shrinkage or /model/stationary: only the "
                   "constant-coefficient Gaussian VARs -- VarNormalWishart, VarNormalGamma and "
                   "VarNormalStochvol -- adapt their prior or restrict their draws that way");
    }
    if (spec.uses_varsel() || spec.structural || spec.n_iid > 0)
    {
        throw std::invalid_argument(
            name + " does not combine /model/shrinkage or /model/stationary with variable "
                   "selection, a structural form or n_iid: each rearranges the prior or the "
                   "draw they act on");
    }
}

/// Refuses a shrinkage prior that does not describe the block: a group per
/// coefficient, groups 1..G with none left out and at least one coefficient
/// shrunk, a diagonal base precision wherever one is, and under `minnesota` a
/// positive shape and rate per group. The starting scales, where given, one per
/// group and one per coefficient, positive.
inline void validate_shrinkage(const Shrinkage type, const ShrinkagePrior &prior,
                               const NormalPrior &base, const arma::vec &initial_scale,
                               const arma::vec &initial_local, const arma::vec &initial_theta,
                               const char *block)
{
    const std::string where = std::string("/priors/") + block + "/shrinkage";
    if (type == Shrinkage::none)
    {
        if (prior.group.n_elem > 0)
        {
            throw std::invalid_argument(where + " is given, and /model/shrinkage is none: name "
                                                "the scheme, or drop the group");
        }
        return;
    }

    const arma::uword n = base.mu.n_elem;
    if (prior.group.n_elem != n)
    {
        throw std::invalid_argument(where + "/group needs one element per coefficient, " +
                                    std::to_string(n) + ", got " +
                                    std::to_string(prior.group.n_elem));
    }
    const arma::uword groups = prior.group.max();
    if (groups == 0)
    {
        throw std::invalid_argument(where + "/group shrinks no coefficient: every element is 0");
    }
    std::vector<bool> seen(groups + 1, false);
    for (arma::uword j = 0; j < n; j++)
    {
        seen[prior.group(j)] = true;
        if (prior.group(j) == 0)
        {
            continue;
        }
        if (!(base.v_inv(j, j) > 0.0))
        {
            throw std::invalid_argument(std::string("/priors/") + block + "/v_inv has a diagonal "
                                        "element that is not positive at shrunk position " +
                                        std::to_string(j + 1));
        }
        for (arma::uword i = 0; i < n; i++)
        {
            if (i != j && (base.v_inv(i, j) != 0.0 || base.v_inv(j, i) != 0.0))
            {
                throw std::invalid_argument(
                    std::string("/priors/") + block + "/v_inv couples shrunk position " +
                    std::to_string(j + 1) + " with position " + std::to_string(i + 1) +
                    "; the scales rescale a diagonal precision only");
            }
        }
    }
    for (arma::uword g = 1; g <= groups; g++)
    {
        if (!seen[g])
        {
            throw std::invalid_argument(where + "/group runs to " + std::to_string(groups) +
                                        " but has no coefficient in group " + std::to_string(g));
        }
    }

    const auto positive = [&](const arma::vec &v, const arma::uword length, const std::string &what) {
        if (v.n_elem != length)
        {
            throw std::invalid_argument(what + " needs " + std::to_string(length) +
                                        " elements, got " + std::to_string(v.n_elem));
        }
        require_finite(v, what);
        if (!(v.min() > 0.0))
        {
            throw std::invalid_argument(what + " must be positive");
        }
    };
    if (type == Shrinkage::minnesota || type == Shrinkage::normal_gamma)
    {
        positive(prior.shape, groups, where + "/shape");
        positive(prior.rate, groups, where + "/rate");
    }
    if (type == Shrinkage::normal_gamma)
    {
        if (prior.theta.n_elem > 0 && prior.theta_rate.n_elem > 0)
        {
            throw std::invalid_argument(where + " has both /theta and /theta_rate: theta is either "
                                                "fixed or drawn, so give one of the two");
        }
        if (prior.theta_rate.n_elem > 0)
        {
            positive(prior.theta_rate, groups, where + "/theta_rate");
            if (initial_theta.n_elem > 0)
            {
                positive(initial_theta, groups, std::string("/initial/") + block + "_theta");
            }
        }
        else
        {
            positive(prior.theta, groups, where + "/theta");
        }
    }
    else if (prior.theta.n_elem > 0 || prior.theta_rate.n_elem > 0)
    {
        throw std::invalid_argument(where + "/theta and /theta_rate are read by /model/shrinkage "
                                            "normal_gamma only, and it is " + to_string(type));
    }
    if (initial_theta.n_elem > 0 && !(type == Shrinkage::normal_gamma && prior.theta_rate.n_elem > 0))
    {
        throw std::invalid_argument(std::string("/initial/") + block + "_theta is read only where "
                                    "/priors/" + block + "/shrinkage/theta_rate draws theta");
    }
    if (initial_scale.n_elem > 0)
    {
        positive(initial_scale, groups, std::string("/initial/") + block + "_shrinkage");
    }
    if (shrinkage_has_local(type) && initial_local.n_elem > 0)
    {
        positive(initial_local, n, std::string("/initial/") + block + "_local");
    }
}

/// Whether the VAR whose k x n_x coefficient matrix is vec'd in `a` is
/// stationary: every eigenvalue of its companion matrix inside the unit
/// circle. The lag block is the first k p columns.
inline bool is_stationary(const arma::vec &a, const arma::uword k, const arma::uword p)
{
    if (p == 0)
    {
        return true;
    }
    const arma::mat a_matrix = arma::reshape(a, k, a.n_elem / k);
    arma::mat companion(k * p, k * p, arma::fill::zeros);
    companion.rows(0, k - 1) = a_matrix.cols(0, k * p - 1);
    if (p > 1)
    {
        companion.submat(k, 0, k * p - 1, k * (p - 1) - 1) = arma::eye<arma::mat>(k * (p - 1), k * (p - 1));
    }
    arma::cx_vec eigenvalues;
    if (!arma::eig_gen(eigenvalues, companion))
    {
        return false;
    }
    return arma::max(arma::abs(eigenvalues)) < 1.0;
}

/// A draw of the coefficients under VarSpec::stationary: `draw` makes one from
/// their conditional posterior, and up to kStationaryTries are made until one
/// is stationary. If none is, the previous draw is kept, which is the reject
/// of an independence Metropolis-Hastings step and so leaves the truncated
/// posterior invariant as the accepted draws do. Returns whether a new draw was
/// kept.
constexpr int kStationaryTries = 100;

inline bool draw_stationary(arma::vec &a, const std::function<arma::vec()> &draw,
                            const arma::uword k, const arma::uword p)
{
    for (int attempt = 0; attempt < kStationaryTries; attempt++)
    {
        arma::vec proposal = draw();
        if (is_stationary(proposal, k, p))
        {
            a = std::move(proposal);
            return true;
        }
    }
    return false;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_SHRINKAGE_SUPPORT_H
