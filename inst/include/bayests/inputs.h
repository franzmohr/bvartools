// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_INPUTS_H
#define BAYESTS_INPUTS_H

#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/spec.h"

namespace bayests
{

/// Where the chain starts.
struct VarNormalWishartInitial
{
    arma::vec a;            ///< Coefficients. Ignored when there are no regressors.
    arma::vec a_lambda;     ///< Inclusion indicators. Ignored without variable selection.
    arma::mat u_sigma_inv;  ///< k x k error precision.
};

/// The complete argument of the VarNormalWishart sampler.
///
/// This is the contract every host fills in: the HDF5 reader builds one out of
/// a file, an R binding builds one out of an S4 object, and neither of them
/// appears anywhere in the numeric code. Adding a field here is how a new
/// input reaches the sampler; there is no other channel.
struct VarNormalWishartInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast; ///< `z` empty when no forecast was requested.

    NormalPrior a_prior;          ///< Unused when there are no regressors.
    WishartPrior u_sigma_prior;
    VarSelPrior varsel_prior;     ///< Unused when spec.varsel is none.

    VarNormalWishartInitial initial;

    /// Whether the model has coefficients to draw at all. Derived from the
    /// data rather than declared, so it cannot disagree with `z`.
    bool use_a() const { return train.nparams() > 0; }

    /// Throws std::invalid_argument describing the first inconsistency it
    /// finds. Called by the sampler before the first draw, so a host that
    /// forgets to call it still gets the message rather than a crash a
    /// thousand iterations in.
    void validate() const;
};

/// Where the VarNormalGamma chain starts.
struct VarNormalGammaInitial
{
    arma::vec a;           ///< Coefficients. Ignored when there are no regressors.
    arma::vec a_lambda;    ///< Inclusion indicators. Ignored without variable selection.
    arma::vec psi;         ///< n_psi free elements of the covariance block.
    arma::vec psi_lambda;
    arma::mat u_sigma_inv; ///< k x k error precision.
};

/// VAR with a normal prior on the coefficients and independent gamma priors on
/// the error precisions, optionally with a constant covariance block.
struct VarNormalGammaInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast; ///< `z` empty when no forecast was requested.

    NormalPrior a_prior;         ///< Unused when there are no regressors.
    VarSelPrior a_varsel_prior;  ///< Unused when spec.varsel is none.

    NormalPrior psi_prior;       ///< Unused without a covariance block.
    VarSelPrior psi_varsel_prior;

    GammaPrior u_sigma_prior;

    VarNormalGammaInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_psi() const { return spec.uses_covar(); }

    void validate() const;
};

/// Where the VarNormalStochvol chain starts.
struct VarNormalStochvolInitial
{
    arma::vec a;
    arma::vec a_lambda;
    arma::vec psi;
    arma::vec psi_lambda;

    arma::mat h;       ///< tt x k log-volatilities, one column per variable.
    arma::vec h_init;  ///< k; the log-volatility before the first period.

    /// k; variance of the log-volatility innovations. A state rather than a
    /// prior -- the sampler redraws it every iteration -- even though the
    /// files keep it next to the prior it is drawn under.
    arma::vec h_sigma;
};

/// VAR with a normal prior on the coefficients and stochastic volatility in
/// the errors, optionally with a covariance block.
///
/// Only BVS reaches this model; SSVS is not implemented for it.
struct VarNormalStochvolInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    NormalPrior a_prior;
    VarSelPrior a_varsel_prior;

    NormalPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    StochvolPrior u_sigma_prior;

    VarNormalStochvolInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_psi() const { return spec.uses_covar(); }

    void validate() const;
};

/// Where the VarTvpGamma chain starts.
///
/// The coefficients are paths rather than points, so most of this is a matrix
/// where the constant-coefficient models hold a vector.
struct VarTvpGammaInitial
{
    arma::mat a;           ///< nparams x tt.
    arma::mat a_sigma_inv; ///< nparams x nparams; only the diagonal is read.
    arma::vec a_init;      ///< nparams; the state before the first period.
    arma::vec a_lambda;

    arma::mat psi;           ///< n_psi x tt.
    arma::mat psi_sigma_inv; ///< n_psi x n_psi; only the diagonal is read.
    arma::vec psi_init;      ///< n_psi.
    arma::vec psi_lambda;

    arma::mat u_omega_inv; ///< k x k.
};

/// VAR whose coefficients follow a random walk, with independent gamma priors
/// on the error precisions and an optional time-varying covariance block.
struct VarTvpGammaInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;

    RandomWalkPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    GammaPrior u_sigma_prior;

    /// Selection for the covariance block. Read from its own group rather than
    /// from the model's, so it can differ from spec.varsel -- and does: this is
    /// the only place where one half of a model selects and the other does not.
    VarSelection psi_varsel = VarSelection::none;

    VarTvpGammaInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_psi() const { return spec.uses_covar(); }
    bool uses_psi_varsel() const { return psi_varsel != VarSelection::none; }

    void validate() const;
};

/// Where the VarTvpWishart chain starts.
///
/// The coefficients are paths rather than points, so most of this is a matrix
/// where the constant-coefficient models hold a vector.
struct VarTvpWishartInitial
{
    arma::mat a;           ///< nparams x tt.
    arma::mat a_sigma_inv; ///< nparams x nparams; only the diagonal is read.
    arma::vec a_init;      ///< nparams; the state before the first period.
    arma::vec a_lambda;

    arma::mat u_sigma_inv;  ///< k x k error precision.
};

/// VAR whose coefficients follow a random walk, with a Wishart prior on the
/// error precision. It carries no psi block: the error covariance is the
/// Wishart precision alone.
struct VarTvpWishartInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;
    WishartPrior u_sigma_prior; 

    VarTvpWishartInitial initial;

    bool use_a() const { return train.nparams() > 0; }

    void validate() const;
};

/// Where the VarTvpStochvol chain starts.
///
/// Both halves of the model carry state: the coefficients are paths, as in
/// VarTvpGamma, and the log-volatility is a path of its own, as in
/// VarNormalStochvol.
struct VarTvpStochvolInitial
{
    arma::mat a;           ///< nparams x tt.
    arma::mat a_sigma_inv; ///< nparams x nparams; only the diagonal is read.
    arma::vec a_init;      ///< nparams; the state before the first period.
    arma::vec a_lambda;

    arma::mat psi;           ///< n_psi x tt.
    arma::mat psi_sigma_inv; ///< n_psi x n_psi; only the diagonal is read.
    arma::vec psi_init;      ///< n_psi.
    arma::vec psi_lambda;

    arma::mat h;      ///< tt x k log-volatilities, one column per variable.
    arma::vec h_init; ///< k; the log-volatility before the first period.

    /// k; variance of the log-volatility innovations. A state rather than a
    /// prior -- the sampler redraws it every iteration -- even though the
    /// files keep it next to the prior it is drawn under.
    arma::vec h_sigma;
};

/// VAR whose coefficients follow a random walk, with stochastic volatility in
/// the errors and an optional time-varying covariance block.
///
/// Only BVS reaches this model; SSVS is not implemented for it.
struct VarTvpStochvolInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;

    RandomWalkPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    StochvolPrior u_sigma_prior;

    /// Selection for the covariance block. Read from its own group rather than
    /// from the model's, so it can differ from spec.varsel, exactly as in
    /// VarTvpGamma.
    VarSelection psi_varsel = VarSelection::none;

    VarTvpStochvolInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_psi() const { return spec.uses_covar(); }
    bool uses_psi_varsel() const { return psi_varsel != VarSelection::none; }

    void validate() const;
};

/// Where the VecNormalWishart chain starts.
struct VecNormalWishartInitial
{
    arma::vec a;
    arma::vec a_lambda;
    arma::vec beta;
    arma::mat u_sigma_inv;
};

/// The complete argument of the VecNormalWishart sampler.
///
/// This is the contract every host fills in: the HDF5 reader builds one out of
/// a file, an R binding builds one out of an S4 object, and neither of them
/// appears anywhere in the numeric code. Adding a field here is how a new
/// input reaches the sampler; there is no other channel.
struct VecNormalWishartInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast; ///< `z` empty when no forecast was requested.

    NormalPrior a_prior;          ///< Unused when there are no regressors.
    VarSelPrior varsel_prior;     ///< Unused when spec.varsel is none.

    ConstantCointSpacePrior beta_prior; ///< Unused when there are no regressors, i.e., r = 0.
    
    WishartPrior u_sigma_prior;    

    VecNormalWishartInitial initial;

    /// Whether the model has coefficients to draw at all.
    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }
    
    /// Throws std::invalid_argument describing the first inconsistency it
    /// finds. Called by the sampler before the first draw, so a host that
    /// forgets to call it still gets the message rather than a crash a
    /// thousand iterations in.
    void validate() const;
};

/// Where the VecNormalGamma chain starts.
struct VecNormalGammaInitial
{
    arma::vec a;
    arma::vec a_lambda;
    arma::vec beta;
    arma::vec psi;         ///< n_psi free elements of the covariance block.
    arma::vec psi_lambda;
    arma::mat u_sigma_inv; ///< k x k error precision.
};

/// VEC with a normal prior on the coefficients, the cointegration space prior on
/// beta and independent gamma priors on the error precisions, optionally with a
/// constant covariance block.
///
/// VecNormalWishart's two coefficient blocks with VarNormalGamma's error block.
/// The one place the swap is visible is Sigma's own posterior: the cointegration
/// space prior conditions alpha on it, and where the Wishart absorbs that back
/// into its scale, independent gammas have no conjugate update for it and
/// bvartools' .bvecalg does not attempt one. Neither does this.
struct VecNormalGammaInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    NormalPrior a_prior;
    VarSelPrior varsel_prior;

    ConstantCointSpacePrior beta_prior; ///< Unused when the rank is zero.

    NormalPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    GammaPrior u_sigma_prior;

    VecNormalGammaInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }
    bool use_psi() const { return spec.uses_covar(); }

    void validate() const;
};

/// Where the VecNormalStochvol chain starts.
struct VecNormalStochvolInitial
{
    arma::vec a;
    arma::vec a_lambda;
    arma::vec beta;
    arma::vec psi;
    arma::vec psi_lambda;

    arma::mat h;      ///< tt x k log-volatilities, one column per variable.
    arma::vec h_init; ///< k; the log-volatility before the first period.
    arma::vec h_sigma;
};

/// VEC with a normal prior on the coefficients, the cointegration space prior on
/// beta and stochastic volatility in the errors, optionally with a constant
/// covariance block.
///
/// Only BVS reaches this model; SSVS is not implemented for it.
///
/// The error precision moves with time and the cointegration space prior needs
/// one to condition alpha on, so the prior takes the average over the sample --
/// the `g_i` of bvartools' .bvecalg. Everything else is VecNormalWishart's
/// coefficient blocks against a per-period precision.
struct VecNormalStochvolInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    NormalPrior a_prior;
    VarSelPrior varsel_prior;

    ConstantCointSpacePrior beta_prior;

    NormalPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    StochvolPrior u_sigma_prior;

    VecNormalStochvolInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }
    bool use_psi() const { return spec.uses_covar(); }

    void validate() const;
};

/// Where the VecTvpWishart chain starts.
///
/// The coefficients and the cointegration vectors are paths; the error
/// precision is not.
struct VecTvpWishartInitial
{
    arma::mat a;           ///< n_a x tt, loadings first.
    arma::mat a_sigma_inv; ///< n_a x n_a; only the diagonal is read.
    arma::vec a_init;      ///< n_a; the state before the first period.
    arma::vec a_lambda;

    arma::mat beta;      ///< n_beta x tt.
    arma::vec beta_init; ///< n_beta.

    arma::mat u_sigma_inv; ///< k x k error precision.
};

/// VEC whose coefficients follow a random walk, cointegration vectors included,
/// with a Wishart prior on the error precision. It carries no psi block: the
/// error covariance is the Wishart precision alone.
struct VecTvpWishartInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;

    TvpCointSpacePrior beta_prior; ///< Unused when the rank is zero.

    WishartPrior u_sigma_prior;

    VecTvpWishartInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }

    void validate() const;
};

/// Where the VecTvpGamma chain starts.
struct VecTvpGammaInitial
{
    arma::mat a;
    arma::mat a_sigma_inv;
    arma::vec a_init;
    arma::vec a_lambda;

    arma::mat beta;
    arma::vec beta_init;

    arma::mat psi;           ///< n_psi x tt.
    arma::mat psi_sigma_inv; ///< n_psi x n_psi; only the diagonal is read.
    arma::vec psi_init;
    arma::vec psi_lambda;

    arma::mat u_omega_inv; ///< k x k.
};

/// VEC whose coefficients follow a random walk, cointegration vectors included,
/// with independent gamma priors on the error precisions and an optional
/// time-varying covariance block.
struct VecTvpGammaInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;

    TvpCointSpacePrior beta_prior;

    RandomWalkPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    GammaPrior u_sigma_prior;

    /// Selection for the covariance block. Read from its own group rather than
    /// from the model's, so it can differ from spec.varsel, exactly as in
    /// VarTvpGamma.
    VarSelection psi_varsel = VarSelection::none;

    VecTvpGammaInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }
    bool use_psi() const { return spec.uses_covar(); }
    bool uses_psi_varsel() const { return psi_varsel != VarSelection::none; }

    void validate() const;
};

/// Where the VecTvpStochvol chain starts.
///
/// Everything this model carries is a path: the coefficients, the cointegration
/// vectors, the covariance block and the log-volatility. `beta` is the one that
/// has no state variance next to it -- see TvpCointSpacePrior for why it is
/// fixed rather than drawn.
struct VecTvpStochvolInitial
{
    arma::mat a;           ///< n_a x tt, loadings first.
    arma::mat a_sigma_inv; ///< n_a x n_a; only the diagonal is read.
    arma::vec a_init;      ///< n_a; the state before the first period.
    arma::vec a_lambda;

    arma::mat beta;      ///< n_beta x tt; each column is vec of a k_beta x rank matrix.
    arma::vec beta_init; ///< n_beta; the state before the first period.

    arma::mat psi;           ///< n_psi x tt.
    arma::mat psi_sigma_inv; ///< n_psi x n_psi; only the diagonal is read.
    arma::vec psi_init;      ///< n_psi.
    arma::vec psi_lambda;

    arma::mat h;      ///< tt x k log-volatilities, one column per variable.
    arma::vec h_init; ///< k; the log-volatility before the first period.

    /// k; variance of the log-volatility innovations. A state rather than a
    /// prior -- the sampler redraws it every iteration -- even though the
    /// files keep it next to the prior it is drawn under.
    arma::vec h_sigma;
};

/// VEC whose loadings, lagged differences and cointegration vectors all follow
/// a random walk, with stochastic volatility in the errors and an optional
/// time-varying covariance block.
///
/// This is bvartools' .bvectvpalg for the "sv" and "sv+covar" error
/// specifications: VarTvpStochvol with one extra Gibbs block, in which beta is
/// drawn as a state path of its own and the loadings' regressors are rebuilt
/// from it. Only BVS reaches this model; SSVS is not implemented for it, and
/// selection may not touch the loadings -- see validate().
struct VecTvpStochvolInput
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast;

    RandomWalkPrior a_prior;
    VarSelPrior a_varsel_prior;

    TvpCointSpacePrior beta_prior; ///< Unused when the rank is zero.

    RandomWalkPrior psi_prior;
    VarSelPrior psi_varsel_prior;

    StochvolPrior u_sigma_prior;

    /// Selection for the covariance block. Read from its own group rather than
    /// from the model's, so it can differ from spec.varsel, exactly as in
    /// VarTvpStochvol.
    VarSelection psi_varsel = VarSelection::none;

    VecTvpStochvolInitial initial;

    bool use_a() const { return train.nparams() > 0; }
    bool use_beta() const { return spec.uses_coint(); }
    bool use_psi() const { return spec.uses_covar(); }
    bool uses_psi_varsel() const { return psi_varsel != VarSelection::none; }

    void validate() const;
};

} // namespace bayests

#endif // BAYESTS_INPUTS_H
