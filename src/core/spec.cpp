// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/spec.h"

#include <limits>
#include <stdexcept>
#include <string>

namespace bayests
{

VarSelection var_selection_from_string(const std::string &name)
{
    if (name == "none" || name.empty())
    {
        return VarSelection::none;
    }
    if (name == "ssvs")
    {
        return VarSelection::ssvs;
    }
    if (name == "bvs")
    {
        return VarSelection::bvs;
    }
    throw std::invalid_argument("unknown variable selection scheme '" + name +
                                "'; expected one of none, ssvs, bvs");
}

const char *to_string(VarSelection selection)
{
    switch (selection)
    {
    case VarSelection::ssvs:
        return "ssvs";
    case VarSelection::bvs:
        return "bvs";
    case VarSelection::none:
        break;
    }
    return "none";
}

Shrinkage shrinkage_from_string(const std::string &name)
{
    if (name == "none" || name.empty())
    {
        return Shrinkage::none;
    }
    if (name == "minnesota")
    {
        return Shrinkage::minnesota;
    }
    if (name == "horseshoe")
    {
        return Shrinkage::horseshoe;
    }
    if (name == "normal_gamma")
    {
        return Shrinkage::normal_gamma;
    }
    throw std::invalid_argument("unknown shrinkage '" + name +
                                "'; expected one of none, minnesota, horseshoe, normal_gamma");
}

const char *to_string(Shrinkage shrinkage)
{
    switch (shrinkage)
    {
    case Shrinkage::minnesota:
        return "minnesota";
    case Shrinkage::horseshoe:
        return "horseshoe";
    case Shrinkage::normal_gamma:
        return "normal_gamma";
    case Shrinkage::none:
        break;
    }
    return "none";
}

ForecastStates forecast_states_from_string(const std::string &name)
{
    if (name == "simulate" || name.empty())
    {
        return ForecastStates::simulate;
    }
    if (name == "hold")
    {
        return ForecastStates::hold;
    }
    throw std::invalid_argument("unknown forecast_states '" + name +
                                "'; expected one of simulate, hold");
}

const char *to_string(ForecastStates states)
{
    return states == ForecastStates::hold ? "hold" : "simulate";
}

UcTrend uc_trend_from_string(const std::string &name)
{
    if (name == "rw" || name.empty())
    {
        return UcTrend::random_walk;
    }
    if (name == "smooth")
    {
        return UcTrend::smooth;
    }
    throw std::invalid_argument("unknown trend '" + name + "'; expected one of rw, smooth");
}

const char *to_string(UcTrend trend)
{
    return trend == UcTrend::smooth ? "smooth" : "rw";
}

UcCycle uc_cycle_from_string(const std::string &name)
{
    if (name == "var" || name.empty())
    {
        return UcCycle::var;
    }
    if (name == "own")
    {
        return UcCycle::own;
    }
    throw std::invalid_argument("unknown cycle '" + name + "'; expected one of var, own");
}

const char *to_string(UcCycle cycle)
{
    return cycle == UcCycle::own ? "own" : "var";
}

UcSkew uc_skew_from_string(const std::string &name)
{
    if (name == "snd" || name.empty())
    {
        return UcSkew::snd;
    }
    if (name == "ald")
    {
        return UcSkew::ald;
    }
    throw std::invalid_argument("unknown skew '" + name + "'; expected one of snd, ald");
}

const char *to_string(UcSkew skew)
{
    return skew == UcSkew::ald ? "ald" : "snd";
}

UcCycleInit uc_cycle_init_from_string(const std::string &name)
{
    if (name == "zero" || name.empty())
    {
        return UcCycleInit::zero;
    }
    if (name == "stationary")
    {
        return UcCycleInit::stationary;
    }
    throw std::invalid_argument("unknown cycle_init '" + name +
                                "'; expected one of zero, stationary");
}

const char *to_string(UcCycleInit init)
{
    return init == UcCycleInit::stationary ? "stationary" : "zero";
}

UcSkewShock uc_skew_shock_from_string(const std::string &name)
{
    if (name == "cycle" || name.empty())
    {
        return UcSkewShock::cycle;
    }
    if (name == "trend")
    {
        return UcSkewShock::trend;
    }
    throw std::invalid_argument("unknown skew_shock '" + name + "'; expected one of cycle, trend");
}

const char *to_string(UcSkewShock shock)
{
    return shock == UcSkewShock::trend ? "trend" : "cycle";
}

void VarSpec::validate() const
{
    if (k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (p < 0 || s < 0 || m < 0 || n < 0 || h < 0)
    {
        throw std::invalid_argument("model dimensions (p, s, m, n, h) cannot be negative");
    }
    if (rank < 0)
    {
        throw std::invalid_argument("cointegration rank cannot be negative");
    }
    if (n_restricted < 0)
    {
        throw std::invalid_argument("the number of restricted deterministic terms cannot be "
                                    "negative");
    }
    if (n_factors < 0 || n_obs_factors < 0)
    {
        throw std::invalid_argument("the number of factors (n_factors, n_obs_factors) cannot be "
                                    "negative");
    }
    if (n_obs_factors > 0 && n_factors == 0)
    {
        throw std::invalid_argument("observed factors (n_obs_factors) belong to a factor augmented "
                                    "VAR, which needs at least one unobserved factor (n_factors) "
                                    "as well; a model with none of the second is a VAR");
    }
    if (rank > k_beta)
    {
        throw std::invalid_argument("cointegration rank (" + std::to_string(rank) +
                                    ") cannot exceed the " + std::to_string(k_beta) +
                                    " rows of the cointegration matrix (k_beta)");
    }
    if (iterations <= 0)
    {
        throw std::invalid_argument("iterations must be positive");
    }
    if (burnin < 0)
    {
        throw std::invalid_argument("burnin cannot be negative");
    }
    if (thin < 1)
    {
        throw std::invalid_argument("thin must be at least 1, which keeps every draw");
    }
    // draws() is an int, and every sampler counts its loop in one.
    if (iterations > (std::numeric_limits<int>::max() - burnin) / thin)
    {
        throw std::invalid_argument(
            "burnin + iterations * thin = " + std::to_string(burnin) + " + " +
            std::to_string(iterations) + " * " + std::to_string(thin) +
            " is longer than a chain can be counted");
    }
}

} // namespace bayests
