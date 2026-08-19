// A refresh of the vendored BayesTS core that loses the Armadillo include
// patch is the one mistake that does not announce itself: upstream's
// `#include <armadillo>` resolves against RcppArmadillo's own copy of the
// Armadillo headers, so everything still compiles, still links and still runs
// -- it just draws from Armadillo's Mersenne twister instead of R's, and
// set.seed() stops reaching the samplers.
//
// This file turns that into a compile error. It includes the core's header
// chain the way the samplers do and asserts that Armadillo came out of it
// configured for R's RNG. See src/core/VENDORED.md.

#include "bayests/inputs.h"

// arma::arma_rng::rng_method is fixed when <armadillo> is first parsed: 2 for
// the alternative RNG that RcppArmadillo installs, 1 or 0 for Armadillo's own.
// Testing ARMA_RNG_ALT instead would not do. One unpatched header is enough to
// configure Armadillo the wrong way, and a patched header later in the same
// translation unit still defines the macro afterwards -- the macro says what
// was asked for, this says what was built.
static_assert(arma::arma_rng::rng_method == 2,
              "The vendored BayesTS core is not using R's RNG. A refresh has "
              "overwritten a \"bayests/arma.h\" include with upstream's "
              "<armadillo>; re-run tools/update-bayests-core.R. "
              "See src/core/VENDORED.md.");
