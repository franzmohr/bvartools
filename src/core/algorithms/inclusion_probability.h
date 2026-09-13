// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_INCLUSION_PROBABILITY_H
#define BAYESTS_CORE_ALGORITHMS_INCLUSION_PROBABILITY_H

#include <cmath>

namespace bayests::core
{

/// The probability of inclusion, exp(l1) / (exp(l0) + exp(l1)), from the log
/// odds l1 - l0 -- the Bernoulli draw both selection schemes end in.
///
/// Both branches are the logistic function; each is the spelling that cannot
/// overflow on its own side. BVS feeds it a log likelihood difference over a
/// whole sample, and SSVS the log ratio of a slab and a spike density, and
/// either is easily in the hundreds. Forming the two weights before dividing,
/// as SSVS used to, sends both to zero at once and the ratio to NaN.
///
/// A prior inclusion probability of exactly zero or one gives log odds of minus
/// or plus infinity, and both branches take those to the limit they should. A
/// NaN comes back as NaN, which a draw `randu() < p` reads as exclusion.
inline double inclusion_probability(const double log_odds)
{
    if (log_odds >= 0.0)
    {
        return 1.0 / (1.0 + std::exp(-log_odds));
    }
    const double odds = std::exp(log_odds);
    return odds / (1.0 + odds);
}

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_INCLUSION_PROBABILITY_H
