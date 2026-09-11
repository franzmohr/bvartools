// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "truncated_normal.h"

#include "bayests/arma.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

/**
 * @file truncated_normal.cpp
 * @brief Draws from a normal distribution restricted to an interval.
 */

namespace bayests::core
{

namespace
{

void require(bool ok, const std::string &what)
{
    if (!ok)
    {
        throw std::invalid_argument("truncated_normal: " + what);
    }
}

/// True if the value is a finite number, read off the exponent bits.
///
/// Not `std::isfinite()`, for the reason `stochvol_mixture.h` and
/// `inverse_gaussian.cpp` set out: a host that compiles these sources with
/// `-ffast-math` is licensed to fold it to `true`, and the guards below are the
/// only thing between a NaN bound and a rejection loop that never ends.
bool is_finite(const double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}

/// Rejection from a uniform over [a, b], with `mode` the point of the interval
/// at which the standard normal density is largest.
///
/// The envelope is that density's value at `mode`, so the acceptance
/// probability is the average of the density over the interval divided by its
/// maximum there. That is close to one for a short interval wherever it sits,
/// which is what makes this the branch for a narrow truncation far out in a
/// tail, where rejection from the untruncated normal would never accept.
double uniform_rejection(const double a, const double b, const double mode)
{
    const double envelope = mode * mode;
    for (;;)
    {
        const double z = a + (b - a) * arma::randu<double>();
        if (arma::randu<double>() <= std::exp(0.5 * (envelope - z * z)))
        {
            return z;
        }
    }
}

/// One draw from the standard normal restricted to [a, b].
double standard_truncated(const double a, const double b)
{
    // Mirror the left tail onto the right one, so the tail branch below only
    // has to handle a lower bound at or above zero. The standard normal is
    // symmetric, so -z is a draw from [-b, -a] whenever z is one from [a, b].
    if (b <= 0.0)
    {
        return -standard_truncated(-b, -a);
    }

    if (a < 0.0)
    {
        // The interval straddles the mode. Rejection from the untruncated
        // normal accepts with probability Phi(b) - Phi(a), which is at least
        // 0.38 once the interval is a standard deviation wide; below that the
        // uniform envelope is the better of the two, and never worse than
        // exp(-1/2).
        if (b - a > 1.0)
        {
            for (;;)
            {
                const double z = arma::randn<double>();
                if (z >= a && z <= b)
                {
                    return z;
                }
            }
        }
        return uniform_rejection(a, b, 0.0);
    }

    // A tail, with the density largest at `a` and falling across the interval.
    // Robert's exponential envelope is shifted to `a` and has effective width
    // 1 / lambda; an interval much shorter than that would spend most of its
    // draws past `b`, so the uniform envelope takes those.
    const double lambda = 0.5 * (a + std::sqrt(a * a + 4.0));
    if ((b - a) * lambda <= 1.0)
    {
        return uniform_rejection(a, b, a);
    }

    for (;;)
    {
        // a + Exp(lambda), by inversion.
        const double z = a - std::log(arma::randu<double>()) / lambda;
        if (z > b)
        {
            continue;
        }
        const double d = z - lambda;
        if (arma::randu<double>() <= std::exp(-0.5 * d * d))
        {
            return z;
        }
    }
}

} // namespace

/**
 * @brief Draws once from a normal distribution restricted to an interval.
 *
 * The draw has density proportional to that of \f$N(\mu, \sigma^2)\f$ on
 * \f$[l, u]\f$ and zero outside it. The interval is standardised to
 * \f$[(l - \mu) / \sigma, (u - \mu) / \sigma]\f$ and one of three rejection
 * schemes draws the standard normal restricted to it:
 *
 * - rejection from the untruncated normal, where the interval contains the mode
 *   and is at least one standard deviation wide;
 * - rejection from a uniform envelope, where the interval is short;
 * - the shifted exponential envelope of Robert (1995), where the interval is a
 *   wide tail.
 *
 * Each of the three accepts with a probability bounded away from zero over the
 * region it is chosen for, so the expected number of iterations is bounded
 * whatever the arguments.
 *
 * @param mu mean of the untruncated distribution, finite.
 * @param sd standard deviation of the untruncated distribution, strictly
 *   positive and finite.
 * @param lower lower end of the interval, finite and below `upper`.
 * @param upper upper end of the interval, finite.
 *
 * @return one draw, in [lower, upper].
 *
 * @throws std::invalid_argument if any argument is not finite, if `sd` is not
 *   strictly positive, or if the interval is empty.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 *
 * Robert, C. P. (1995). Simulation of truncated normal variables. Statistics
 * and Computing, 5(2), 121-125.
 */
double truncated_normal(const double mu, const double sd, const double lower, const double upper)
{
    require(is_finite(mu), "'mu' is NaN or infinite");
    require(is_finite(sd), "'sd' is NaN or infinite");
    require(is_finite(lower), "'lower' is NaN or infinite");
    require(is_finite(upper), "'upper' is NaN or infinite");
    require(sd > 0.0, "'sd' must be strictly positive");
    require(lower < upper, "'lower' must be below 'upper'");

    const double z = standard_truncated((lower - mu) / sd, (upper - mu) / sd);

    // Clamped rather than returned as it stands: the draw is inside the
    // standardised interval by construction, but undoing the standardisation is
    // two rounded operations and can land a hair outside a bound that a caller
    // is about to treat as a hard constraint.
    const double out = mu + sd * z;
    return out < lower ? lower : (out > upper ? upper : out);
}

} // namespace bayests::core
