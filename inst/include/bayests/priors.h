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

/// Uniform prior on the autoregression of a cointegration space that moves with
/// time, and the switch that turns the draw of it on.
///
/// Koop, Leon-Gonzalez and Strachan (2011) treat rho as a parameter and put a
/// uniform prior on a range close to one -- (0.999, 1) in their application --
/// on the grounds that it is rho, the innovation variance being fixed at the
/// identity, that controls how concentrated the distribution of the
/// cointegration space at t is around the space at t - 1. The support is left
/// to the file rather than fixed here, but the default is theirs.
///
/// Off unless a file asks for it, so a model file that does not name the
/// support runs with the fixed `TvpCointSpacePrior::rho` it always did.
struct CointRhoPrior
{
    /// Whether rho is drawn. The two bounds below mean nothing when it is not.
    bool draw = false;

    double min = 0.999; ///< Lower end of the support, above zero.
    double max = 1.0;   ///< Upper end, at most one.
};

/// Prior on a cointegration space that moves with time:
///
///     beta_t = rho (I_r kron P_tau) beta_{t-1} + eta_t,   eta_t ~ N(0, I).
///
/// With `p_tau` empty the transition is rho alone: the space at t is centred on
/// the space at t - 1, and its marginal prior is uniform at every t -- Koop,
/// Leon-Gonzalez and Strachan (2011, eq. 6). With `p_tau` set it is their
/// informative marginal prior instead (working paper version, eq. 12).
///
/// The innovation variance is the identity and is deliberately not a knob. Only
/// the product alpha beta' is identified, so something has to fix beta's scale;
/// where the constant-coefficient VEC normalises the draw after the fact, this
/// one lets the state equation do it -- the same normalisation bvartools'
/// .bvectvpalg makes by hardcoding a unit state variance. The scale of the
/// relation then lives in alpha, whose own state variance is drawn.
struct TvpCointSpacePrior
{
    /// Autoregression of the state equation: the value the chain starts at, and
    /// the value it keeps for the whole run unless `rho_prior` turns its draw
    /// on.
    ///
    /// The default is Koop, Leon-Gonzalez and Strachan's rather than the random
    /// walk `.bvectvpalg` hardcodes. Just below one, beta_t has the stationary
    /// distribution N(0, I / (1 - rho^2)), so the prior is proper and the path
    /// is pulled back towards the space `initial_state` names; at exactly one it
    /// is a random walk, whose variance grows without bound over the sample and
    /// which, beta being identified only up to scale, has nothing to pull it
    /// back. One is still accepted, and is what a file that predates this field
    /// means.
    double rho = 0.999;

    /// Uniform prior on `rho`, and the switch that turns its draw on.
    CointRhoPrior rho_prior;

    /// k_beta x k_beta: the transition of the state equation with rho taken out.
    ///
    /// Koop, Leon-Gonzalez and Strachan's informative marginal prior sets it to
    /// P_tau = H H' + tau H_perp H_perp' with H semi-orthogonal and 0 <= tau <= 1.
    /// The part of beta along sp(H) then keeps rho and the part off it decays at
    /// rho tau, so the space at t is centred between the space at t - 1 and sp(H),
    /// and the mode of its marginal distribution is sp(H) at every t. It enters as
    /// the transition rather than the innovation variance, which is what leaves
    /// beta's scale pinned as above. Any symmetric matrix with eigenvalues in
    /// [0, 1] is accepted -- a tau per direction rather than one.
    ///
    /// Empty means the identity: the noninformative prior, and what every file
    /// written before this field describes.
    ///
    /// Their prior on the state at the start of the sample is the stationary
    /// distribution the transition implies, N(0, I_r kron P_tau* / (1 - rho^2))
    /// with tau* = (1 - rho^2) / (1 - rho^2 tau^2). Here that is `initial_state`
    /// below, which the file supplies -- for this prior as for the identity.
    arma::mat p_tau;

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
