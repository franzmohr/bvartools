// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_ARMA_H
#define BAYESTS_ARMA_H

/**
 * @file arma.h
 * @brief The single point at which Armadillo enters BayesTS.
 *
 * Every header and translation unit that needs Armadillo includes this rather
 * than `<armadillo>` directly, so that a host which has to supply its own
 * Armadillo can redirect all of them at once instead of patching each file.
 *
 * A standalone build needs nothing: this defaults to `<armadillo>`. An embedded
 * host that must be seen first defines the macro instead --
 *
 *     target_compile_definitions(... PRIVATE
 *         BAYESTS_ARMA_HEADER=<RcppArmadillo.h>)
 *
 * or, for a build system that is not CMake,
 *
 *     -DBAYESTS_ARMA_HEADER='<RcppArmadillo.h>'
 *
 * The R case is the reason this exists. RcppArmadillo has to be the first thing
 * every translation unit sees, because it is what points Armadillo's RNG at R's
 * own: include `<armadillo>` first and the samplers draw from Armadillo's
 * Mersenne twister instead, which `set.seed()` cannot reach and `R CMD check`
 * does not permit.
 *
 * @warning A host that vendors these sources and does *not* define
 * `BAYESTS_ARMA_HEADER` gets the standalone default, which compiles and runs
 * but silently draws from the wrong random number generator. The failure is not
 * a build error, so it is worth asserting in the host that the macro is set.
 */

#ifndef BAYESTS_ARMA_HEADER
#define BAYESTS_ARMA_HEADER <armadillo>
#endif

#include BAYESTS_ARMA_HEADER

#endif // BAYESTS_ARMA_H
