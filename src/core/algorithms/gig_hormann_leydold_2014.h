// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_GIG_HORMANN_LEYDOLD_2014_H
#define BAYESTS_CORE_ALGORITHMS_GIG_HORMANN_LEYDOLD_2014_H

namespace bayests::core
{

/// One draw from the generalised inverse Gaussian distribution GIG(lambda,
/// chi, psi), with density proportional to
///
///     x^(lambda - 1) exp(-(chi / x + psi x) / 2),   x > 0.
///
/// The three samplers of Hoermann and Leydold (2014), each where it has a
/// uniformly bounded rejection constant: the ratio-of-uniforms method around the
/// mode for lambda > 2 or sqrt(chi psi) > 3, the one without the shift where
/// lambda >= 1 - 2.25 chi psi or sqrt(chi psi) > 0.2, and their rejection
/// method for the density that is not T-concave everywhere else. So no
/// combination of parameters makes it slow, which matters to a caller whose
/// parameters are those of the last draw.
///
/// chi = 0 is the gamma limit and needs lambda > 0; psi = 0 the inverse gamma
/// one and needs lambda < 0. Both zero is refused, as is a negative or
/// non-finite argument.
///
/// The number of variates consumed per draw is not fixed. A chain that reaches
/// this is still reproducible under a fixed seed, but not comparable with one
/// whose arguments differed.
///
/// Hoermann, W., & Leydold, J. (2014). Generating generalized inverse Gaussian
/// random variates. Statistics and Computing, 24(4), 547-557.
double gig_hormann_leydold_2014(double lambda, double chi, double psi);

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_GIG_HORMANN_LEYDOLD_2014_H
