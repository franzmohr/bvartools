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

/// Prior on a cointegration space that moves with time:
///
///     beta_t = rho beta_{t-1} + eta_t,   eta_t ~ N(0, I).
///
/// The innovation variance is the identity and is deliberately not a knob. Only
/// the product alpha beta' is identified, so something has to fix beta's scale;
/// where the constant-coefficient VEC normalises the draw after the fact, this
/// one lets the state equation do it -- the same normalisation bvartools'
/// .bvectvpalg makes by hardcoding a unit state variance. The scale of the
/// relation then lives in alpha, whose own state variance is drawn.
struct TvpCointSpacePrior
{
    /// Autoregression of the state equation, and not drawn -- Koop,
    /// Leon-Gonzalez and Strachan (2011) sample it, bvartools does not, and
    /// neither does this.
    ///
    /// The default is theirs rather than the random walk `.bvectvpalg`
    /// hardcodes. Just below one, beta_t has the stationary distribution
    /// N(0, I / (1 - rho^2)), so the prior is proper and the path is pulled back
    /// towards the space `initial_state` names; at exactly one it is a random
    /// walk, whose variance grows without bound over the sample and which,
    /// beta being identified only up to scale, has nothing to pull it back.
    /// One is still accepted, and is what a file that predates this field means.
    double rho = 0.999;

    /// Normal on the state of the period before the sample.
    NormalPrior initial_state;
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
