// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/spec.h"

#include <stdexcept>

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
}

} // namespace bayests
