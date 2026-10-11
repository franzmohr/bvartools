// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_ERROR_SCALE_CHECKS_H
#define BAYESTS_CORE_MODELS_ERROR_SCALE_CHECKS_H

// The refusals of VarNormalWishart's two common scales on the error
// covariance, Student t errors and a volatility break, without the samplers
// that draw them (student_t_support.h, volatility_break_support.h).
//
// Every algorithm's validate() makes the require_supported_*() calls and
// VarNormalWishart's makes the validate_*() ones, so src/core/inputs.cpp needs
// these and nothing else of either model. Kept apart so that a package
// vendoring the core by its include closure -- a factor or unobserved
// components subset -- does not pull in a sampler block none of its models
// runs.

#include "bayests/arma.h"
#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/spec.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Refuses /model/student_t where the algorithm does not read it.
inline void require_supported_student_t(const VarSpec &spec, const bool supported,
                                        const char *algorithm)
{
    if (spec.student_t && !supported)
    {
        throw std::invalid_argument(std::string(algorithm) +
                                    " does not read /model/student_t: only VarNormalWishart has "
                                    "Student t errors");
    }
}

/// Refuses Student t errors the sampler cannot run as given: a prior that is
/// neither a fixed nu nor a gamma on it, or both; a start that is not
/// positive; the steady-state prior, whose draw of the mean reads the error
/// covariance as the same in every period; a score from a panel or of a horizon
/// not observed whole, where lambda does not integrate out in closed form; and
/// a scenario, which carries information about the horizon's lambda that a
/// draw from its prior ignores. A training panel not observed whole is read:
/// its completion takes one covariance per period, and its log likelihood
/// conditions on the lambdas drawn.
inline void validate_student_t(const VarSpec &spec, const StudentTPrior &prior,
                               const double initial_df, const TrainData &train,
                               const ForecastData &forecast, const TestData &test,
                               const char *algorithm)
{
    if (!spec.student_t)
    {
        return;
    }
    const std::string name(algorithm);
    const bool fixed = prior.df != 0.0;
    const bool drawn = prior.shape != 0.0 || prior.rate != 0.0;
    if (fixed == drawn)
    {
        throw std::invalid_argument(
            name + ": /model/student_t needs exactly one of /priors/student_t/df, which holds nu "
                   "fixed, and /priors/student_t/shape and rate, a gamma prior that draws it");
    }
    const auto positive = [&](const double v, const char *what) {
        if (!(v > 0.0) || !std::isfinite(v))
        {
            throw std::invalid_argument(name + ": /priors/student_t/" + what +
                                        " must be positive and finite; got " + std::to_string(v));
        }
    };
    if (fixed)
    {
        positive(prior.df, "df");
    }
    else
    {
        positive(prior.shape, "shape");
        positive(prior.rate, "rate");
    }
    if (initial_df != 0.0 && (!(initial_df > 0.0) || !std::isfinite(initial_df)))
    {
        throw std::invalid_argument(name + ": /initial/student_t_df must be positive and finite; got " +
                                    std::to_string(initial_df));
    }
    if (spec.steady_state)
    {
        throw std::invalid_argument(
            name + " does not combine /model/student_t with /model/steady_state, whose draw of "
                   "the unconditional mean reads the error covariance as the same in every period");
    }
    if (!test.constraints.empty() || (!train.constraints.empty() && test.y.n_elem > 0))
    {
        throw std::invalid_argument(
            name + " does not score a forecast with Student t errors from a panel not observed "
                   "whole, or against a horizon realised in part (/data/test/constraints), yet: the "
                   "horizon's lambda would have to be integrated out of a density whose rows reach "
                   "across periods");
    }
    if (!forecast.constraints.empty())
    {
        throw std::invalid_argument(
            name + " does not condition a forecast with Student t errors on "
                   "/data/forecast/constraints yet: a scenario says something about how large the "
                   "horizon's shocks are, so lambda would have to be drawn given it, not from its "
                   "prior");
    }
}

/// Refuses /model/volatility_break where the algorithm does not read it.
inline void require_supported_volatility_break(const VarSpec &spec, const bool supported,
                                               const char *algorithm)
{
    if (spec.volatility_break != 0 && !supported)
    {
        throw std::invalid_argument(
            std::string(algorithm) +
            " does not read /model/volatility_break: only VarNormalWishart scales its error "
            "covariance from a known period on");
    }
}

/// Refuses a volatility break the sampler cannot estimate as given: a break
/// outside the sample, no free scale, a prior that is not proper, a start
/// outside the support, and the steady-state prior, whose draw of the mean
/// reads the error covariance as the same in every period. A panel or a
/// horizon not observed whole is read: the completion step takes one
/// covariance per period.
inline void validate_volatility_break(const VarSpec &spec, const VolatilityBreakPrior &prior,
                                      const arma::vec &initial, const arma::uword tt,
                                      const char *algorithm)
{
    const std::string name(algorithm);
    if (spec.volatility_break < 0)
    {
        throw std::invalid_argument(name + ": /model/volatility_break is a period of "
                                           "/data/train/y counted from one, or zero for none; got " +
                                    std::to_string(spec.volatility_break));
    }
    if (!spec.uses_volatility_break())
    {
        return;
    }
    if (spec.volatility_scales < 1)
    {
        throw std::invalid_argument(name + ": /model/volatility_scales must be at least one, the "
                                           "period of the break itself; got " +
                                    std::to_string(spec.volatility_scales));
    }
    if (static_cast<arma::uword>(spec.volatility_break) > tt)
    {
        throw std::invalid_argument(
            name + ": /model/volatility_break is period " + std::to_string(spec.volatility_break) +
            " and /data/train/y has " + std::to_string(tt) +
            ": the break has to be inside the sample, whose periods after it are what "
            "estimate it");
    }
    if (spec.steady_state)
    {
        throw std::invalid_argument(
            name + " does not combine /model/volatility_break with /model/steady_state, whose "
                   "draw of the unconditional mean reads the error covariance as the same in "
                   "every period");
    }
    const auto positive = [&](const double v, const char *what) {
        if (!(v > 0.0) || !std::isfinite(v))
        {
            throw std::invalid_argument(name + ": /priors/volatility/" + what +
                                        " must be positive and finite under "
                                        "/model/volatility_break; got " +
                                        std::to_string(v));
        }
    };
    positive(prior.shape, "shape");
    positive(prior.rho_shape1, "rho_shape1");
    positive(prior.rho_shape2, "rho_shape2");

    if (initial.n_elem > 0)
    {
        const arma::uword m = static_cast<arma::uword>(spec.volatility_scales);
        if (initial.n_elem != m + 1)
        {
            throw std::invalid_argument(name + ": /initial/volatility needs volatility_scales + 1 = " +
                                        std::to_string(m + 1) +
                                        " elements, the free scales and then the decay rate; got " +
                                        std::to_string(initial.n_elem));
        }
        for (arma::uword j = 0; j < m; j++)
        {
            if (!(initial(j) > 1.0) || !std::isfinite(initial(j)))
            {
                throw std::invalid_argument(name + ": every scale in /initial/volatility must be "
                                                   "finite and above one; element " +
                                            std::to_string(j + 1) + " is " +
                                            std::to_string(initial(j)));
            }
        }
        if (!(initial(m) > 0.0 && initial(m) < 1.0))
        {
            throw std::invalid_argument(name + ": the decay rate, the last element of "
                                               "/initial/volatility, must be in (0, 1); got " +
                                        std::to_string(initial(m)));
        }
    }
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_ERROR_SCALE_CHECKS_H
