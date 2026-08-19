// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_PRIORS_H
#define BAYESTS_PRIORS_H

#include "bayests/arma.h"

namespace bayests
{

/// Normal prior on a coefficient vector.
///
/// The precision is stored rather than the covariance because that is what the
/// posterior update needs, and because SSVS rewrites its diagonal in place
/// once per draw.
struct NormalPrior
{
    arma::vec mu;
    arma::mat v_inv;
};

/// Prior on the cointegration space.
struct ConstantCointSpacePrior
{
    double v_inv;
    arma::mat p_tau_inv;
};

/// Wishart prior on an inverse covariance matrix.
struct WishartPrior
{
    int df = 0;
    arma::mat scale;
};

/// Independent gamma priors, one per element.
///
/// Rate rather than scale, because that is what the files store and what the
/// posterior update adds the sum of squares to. Armadillo's randg() wants a
/// scale, so the samplers invert at the point of use.
struct GammaPrior
{
    arma::vec shape;
    arma::vec rate;
};

/// A random walk state equation: how far the state may drift from one period
/// to the next, and where it starts.
struct RandomWalkPrior
{
    /// Inverse gamma on the variance of the state innovations.
    GammaPrior sigma;

    /// Normal on the state of the period before the sample.
    NormalPrior initial_state;
};

/// Everything the stochastic volatility block reads beyond the state equation.
struct StochvolPrior
{
    /// Added inside the log before the mixture approximation is applied, so a
    /// residual that lands on zero does not send log(u^2) to -inf and take the
    /// whole chain with it.
    arma::vec offset;

    /// The log-volatility follows a random walk of its own.
    RandomWalkPrior state;
};

/// The half of a variable selection prior that SSVS needs: the two component
/// standard deviations of the spike-and-slab mixture.
struct SsvsPrior
{
    arma::vec tau0; ///< Spike; the "excluded" component.
    arma::vec tau1; ///< Slab; the "included" component.
};

/// Everything the selection step reads, for either scheme.
struct VarSelPrior
{
    /// Prior inclusion probability, one per coefficient.
    arma::vec inprior;

    /// Zero-based positions of the coefficients selection applies to. The
    /// HDF5 files and R both count from one; conversion happens at the edge so
    /// the sampler never has to remember which convention it is holding.
    arma::uvec include;

    SsvsPrior ssvs;
    /// Not needed for BVS, because it only needs the positions of the coefficients
    /// that are selected, and those are already in `include`.
    /// The sampler does not need to know which coefficients are excluded, because it
    /// only needs to know which ones are included.

    arma::uword size() const { return include.n_elem; }
};

} // namespace bayests

#endif // BAYESTS_PRIORS_H
