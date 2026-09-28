// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "gig_hormann_leydold_2014.h"
#include "bayests/arma.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

/**
 * @file gig_hormann_leydold_2014.cpp
 * @brief Draws from the generalised inverse Gaussian distribution.
 *
 * Everything below works with the standardised density
 *
 *     f(x) = x^(l - 1) exp(-w (x + 1/x) / 2),   l >= 0, w > 0,
 *
 * which GIG(lambda, chi, psi) is a rescaling of: with w = sqrt(chi psi) and
 * alpha = sqrt(chi / psi), alpha X is GIG(|lambda|, chi, psi) for X drawn from
 * f with l = |lambda|, and alpha / X is GIG(-|lambda|, chi, psi). Densities are
 * handled as logarithms relative to their value at the mode, so that neither a
 * tiny w -- a mode near zero and a tail out at 2 / w -- nor a large one
 * overflows.
 */

namespace bayests::core
{

namespace
{

/// True if x is a finite number, read off the exponent bits: not
/// std::isfinite(), which a host compiling with -ffast-math may fold to true.
bool is_finite(const double x)
{
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}

/// A uniform on (0, 1), never zero, so that its logarithm is finite.
double positive_uniform()
{
    double u = 0.0;
    while (u <= 0.0)
    {
        u = arma::randu<double>();
    }
    return u;
}

/// The mode of x^(l - 1) exp(-w (x + 1/x) / 2): the positive root of
/// w x^2 / 2 - (l - 1) x - w / 2. Written as w / (r - a) where a < 0, which is
/// the same number without the cancellation of (a + r) / w.
double mode(const double l, const double w)
{
    const double a = l - 1.0;
    const double r = std::hypot(a, w);
    return a >= 0.0 ? (a + r) / w : w / (r - a);
}

/// log f(x) - log f(m).
double log_ratio(const double x, const double m, const double l, const double w)
{
    return (l - 1.0) * std::log(x / m) - 0.5 * w * (x + 1.0 / x - m - 1.0 / m);
}

/// Ratio of uniforms without a shift, for the region where Hoermann and
/// Leydold show its rejection constant is bounded. v runs to sqrt(f(m)) and u
/// to the maximum of x sqrt(f(x)), which is the mode of x^(l + 1) exp(...).
double ratio_of_uniforms(const double l, const double w)
{
    const double m = mode(l, w);
    const double x_u = mode(l + 2.0, w);
    const double u_max = x_u * std::exp(0.5 * log_ratio(x_u, m, l, w));
    for (;;)
    {
        const double u = u_max * arma::randu<double>();
        const double v = positive_uniform();
        const double x = u / v;
        if (x > 0.0 && 2.0 * std::log(v) <= log_ratio(x, m, l, w))
        {
            return x;
        }
    }
}

/// Ratio of uniforms around the mode, for l > 2 or w > 3. The bounds on u are
/// the extremes of (x - m) sqrt(f(x)) on either side of the mode, where the
/// derivative of its logarithm,
///
///     1 / (x - m) + ((l - 1) / x - w / 2 + w / (2 x^2)) / 2,
///
/// changes sign from positive to negative. It is monotone in between, so
/// bisection finds both to the last digit without the cubic's closed form and
/// its cancellations.
double shifted_ratio_of_uniforms(const double l, const double w)
{
    const double m = mode(l, w);
    const auto slope = [&](const double x) {
        return 1.0 / (x - m) + 0.5 * ((l - 1.0) / x - 0.5 * w + 0.5 * w / (x * x));
    };
    const auto bisect = [&](double lo, double hi) {
        for (int i = 0; i < 200 && hi - lo > 4.0 * 2.220446049250313e-16 * hi; i++)
        {
            const double mid = 0.5 * (lo + hi);
            (slope(mid) > 0.0 ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };

    const double x_lo = bisect(0.0, m);
    double hi = 2.0 * m + 1.0;
    while (slope(hi) > 0.0)
    {
        hi *= 2.0;
    }
    const double x_hi = bisect(m, hi);

    const double u_min = (x_lo - m) * std::exp(0.5 * log_ratio(x_lo, m, l, w));
    const double u_max = (x_hi - m) * std::exp(0.5 * log_ratio(x_hi, m, l, w));
    for (;;)
    {
        const double u = u_min + (u_max - u_min) * arma::randu<double>();
        const double v = positive_uniform();
        const double x = u / v + m;
        if (x > 0.0 && 2.0 * std::log(v) <= log_ratio(x, m, l, w))
        {
            return x;
        }
    }
}

/// The rejection method of Hoermann and Leydold for l < 1 and small w, where
/// the density is not T-concave. The hat is f(m) on (0, x0), exp(-w) x^(l - 1)
/// on (x0, xs), since x + 1/x >= 2, and xs^(l - 1) exp(-w x / 2) past xs, since
/// x^(l - 1) falls and exp(-w / (2 x)) is below one; x0 = w / (1 - l) and
/// xs = max(x0, 2 / w). Each piece is drawn by inversion.
double non_t_concave(const double l, const double w)
{
    const double m = mode(l, w);
    const double log_fm = (l - 1.0) * std::log(m) - 0.5 * w * (m + 1.0 / m);
    const double x0 = w / (1.0 - l);
    const double xs = std::max(x0, 2.0 / w);
    const double log_x0 = std::log(x0);
    const double log_xs = std::log(xs);
    const double span = log_xs - log_x0;

    // The area under each piece, relative to f(m) and as logarithms.
    const double log_a1 = log_x0;
    double log_a2 = -INFINITY;
    if (span > 0.0)
    {
        // The integral of x^(l - 1) over (x0, xs), x0^l expm1(l span) / l, and
        // span itself in the limit l = 0.
        const double integral =
            l > 0.0 ? l * log_x0 + std::log(std::expm1(l * span) / l) : std::log(span);
        log_a2 = -w - log_fm + integral;
    }
    const double log_a3 = (l - 1.0) * log_xs - log_fm + std::log(2.0 / w) - 0.5 * w * xs;

    const double top = std::max({log_a1, log_a2, log_a3});
    const double a1 = std::exp(log_a1 - top);
    const double a2 = std::exp(log_a2 - top);
    const double a3 = std::exp(log_a3 - top);
    const double total = a1 + a2 + a3;

    for (;;)
    {
        const double pick = total * arma::randu<double>();
        const double log_v = std::log(positive_uniform());
        if (pick < a1)
        {
            const double x = x0 * positive_uniform();
            if (log_v <= (l - 1.0) * std::log(x) - 0.5 * w * (x + 1.0 / x) - log_fm)
            {
                return x;
            }
        }
        else if (pick < a1 + a2)
        {
            const double u = arma::randu<double>();
            const double log_x =
                l > 0.0 ? log_x0 + std::log1p(u * std::expm1(l * span)) / l : log_x0 + u * span;
            const double x = std::exp(log_x);
            if (log_v <= w - 0.5 * w * (x + 1.0 / x))
            {
                return x;
            }
        }
        else
        {
            const double x = xs - 2.0 / w * std::log(positive_uniform());
            if (log_v <= (l - 1.0) * (std::log(x) - log_xs) - 0.5 * w / x)
            {
                return x;
            }
        }
    }
}

/// One draw from the standardised density with l >= 0 and w > 0.
double standard_gig(const double l, const double w)
{
    if (l > 2.0 || w > 3.0)
    {
        return shifted_ratio_of_uniforms(l, w);
    }
    if (l >= 1.0 - 2.25 * w * w || w > 0.2)
    {
        return ratio_of_uniforms(l, w);
    }
    return non_t_concave(l, w);
}

} // namespace

double gig_hormann_leydold_2014(const double lambda, const double chi, const double psi)
{
    if (!is_finite(lambda) || !is_finite(chi) || !is_finite(psi) || !(chi >= 0.0) ||
        !(psi >= 0.0))
    {
        throw std::invalid_argument("gig_hormann_leydold_2014: needs a finite lambda and finite, "
                                    "non-negative chi and psi; got lambda = " +
                                    std::to_string(lambda) + ", chi = " + std::to_string(chi) +
                                    ", psi = " + std::to_string(psi));
    }

    const double w = std::sqrt(chi) * std::sqrt(psi);
    if (w > 0.0)
    {
        const double alpha = std::sqrt(chi) / std::sqrt(psi);
        const double x = standard_gig(std::fabs(lambda), w);
        return lambda < 0.0 ? alpha / x : alpha * x;
    }

    // The limits: a gamma with rate psi / 2 where chi vanishes, the reciprocal
    // of one with rate chi / 2 where psi does. Each is proper only on its own
    // side of lambda = 0.
    if (lambda > 0.0 && psi > 0.0)
    {
        return arma::randg<double>(arma::distr_param(lambda, 2.0 / psi));
    }
    if (lambda < 0.0 && chi > 0.0)
    {
        return 1.0 / arma::randg<double>(arma::distr_param(-lambda, 2.0 / chi));
    }
    throw std::invalid_argument("gig_hormann_leydold_2014: GIG(" + std::to_string(lambda) + ", " +
                                std::to_string(chi) + ", " + std::to_string(psi) +
                                ") is improper: chi = 0 needs lambda > 0 and psi > 0, psi = 0 "
                                "needs lambda < 0 and chi > 0");
}

} // namespace bayests::core
