// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_TRUNCATED_NORMAL_H
#define BAYESTS_CORE_ALGORITHMS_TRUNCATED_NORMAL_H

namespace bayests::core
{

/// One draw from N(mu, sd^2) restricted to the interval [lower, upper].
///
/// Exact -- rejection sampling from the untruncated law and from two envelopes
/// that cover it, never an approximation of the truncated quantile function,
/// which loses all its digits once the interval sits far out in a tail. That
/// case is not hypothetical here: the caller is the autoregression of a
/// cointegration state equation, whose prior support is a narrow interval below
/// one and whose conditional mean can land far outside it, so the standardised
/// interval is routinely tens of standard deviations from the mode.
///
/// The number of variates consumed per draw is therefore not fixed. A chain
/// that reaches this is still reproducible under a fixed seed -- the sequence
/// is a function of the seed and the arguments -- but it is not comparable with
/// one whose arguments differed.
double truncated_normal(double mu, double sd, double lower, double upper);

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_TRUNCATED_NORMAL_H
