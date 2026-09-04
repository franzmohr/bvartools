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

/// Where the DfmNormalGamma chain starts.
///
/// The factor path is not among these: it is the first thing every draw
/// produces, from the loadings and the two precisions, so there is nothing for a
/// starting value to do.
struct DfmNormalGammaInitial
{
    /// The `n_lambda` free loadings, in the order the sampler draws them --
    /// row by row, and within a row left to right. See DfmNormalGammaInput.
    arma::vec lambda;

    /// vec([A_1 .. A_p]). Ignored when the factors have no dynamics.
    arma::vec a;

    /// k; the diagonal of the idiosyncratic precision. A vector rather than a
    /// matrix because U is diagonal by assumption and nothing here would read
    /// an off-diagonal element.
    arma::vec u_sigma_inv;

    /// n_factors; the diagonal of the factor innovation precision.
    arma::vec v_sigma_inv;
};

/// Dynamic factor model with a normal prior on the loadings and on the factor
/// transition, and independent gamma priors on both error precisions.
///
///     x_t = Lambda f_t + u_t,                  u_t ~ N(0, U),  U diagonal,
///     f_t = sum_{j=1..p} A_j f_{t-j} + v_t,    v_t ~ N(0, V),  V diagonal.
///
/// `train.y` holds the observed series, tt x k, and is the only data the model
/// takes -- a DFM has no regressors, so `z`, `x` and `w` are all unused. There is
/// no ForecastData either: the forecast is a simulation of the transition
/// forward from the last drawn factors, which needs no out-of-sample matrix.
/// `spec.h` is the horizon.
///
/// Two orderings have to be agreed with the host and are the only ones that are
/// not implied by a dimension:
///
///   - the free loadings, everywhere they appear as a vector -- `initial.lambda`
///     and both halves of `lambda_prior` -- run row by row: row 1's single free
///     element, then row 2's two, and so on, each row left to right. That is the
///     order the equation-by-equation draw consumes them in, which is what makes
///     a row's prior block contiguous. Note that the *posterior* is stored
///     differently, as vec of the whole M x N matrix; see DfmNormalGammaDraws.
///   - `a` is vec([A_1 .. A_p]) with the blocks side by side, so the first
///     n_factors^2 elements are A_1 read column-wise.
///
/// Values of the factors before the sample are zero rather than drawn, which is
/// what makes the first p periods' transitions well defined and is the
/// convention bvartools' dfmpost() uses.
struct DfmNormalGammaInput
{
    VarSpec spec;
    TrainData train;

    NormalPrior lambda_prior;  ///< Over the free loadings, in the row-major order above.
    NormalPrior a_prior;       ///< Unused when the factors have no dynamics.

    GammaPrior u_sigma_prior;  ///< k independent priors on the idiosyncratic precisions.
    GammaPrior v_sigma_prior;  ///< n_factors independent priors on the factor innovations.

    DfmNormalGammaInitial initial;

    /// Whether the factors carry dynamics to draw. A DFM of transition order
    /// zero is a static factor model with serially independent factors, which
    /// this sampler estimates -- it is not an error.
    bool use_a() const { return spec.n_factor_a() > 0; }

    /// Throws std::invalid_argument describing the first inconsistency it
    /// finds. Called by the sampler before the first draw, so a host that
    /// forgets to call it still gets the message rather than a crash a
    /// thousand iterations in.
    void validate() const;
};

/// Where the DfmNormalStochvol chain starts.
///
/// DfmNormalGammaInitial plus a log-volatility path for each of the two error
/// terms, and minus the two precisions those paths replace. As there, the factor
/// path itself is not here: it is the first thing every draw produces.
struct DfmNormalStochvolInitial
{
    /// The `n_lambda` free loadings, row by row and within a row left to right.
    /// See DfmNormalStochvolInput.
    arma::vec lambda;

    /// vec([A_1 .. A_p]). Ignored when the factors have no dynamics.
    arma::vec a;

    /// tt x k log-volatilities of the idiosyncratic errors, one column per
    /// observed series.
    arma::mat u_h;

    /// k; the idiosyncratic log-volatility before the first period.
    arma::vec u_h_init;

    /// k; variance of the idiosyncratic log-volatility innovations. A state
    /// rather than a prior -- the sampler redraws it every iteration -- even
    /// though the files keep it next to the prior it is drawn under, which is
    /// the convention VarNormalStochvolInitial already set.
    arma::vec u_h_sigma;

    /// tt x n_factors log-volatilities of the factor innovations.
    arma::mat v_h;

    /// n_factors; the factor-innovation log-volatility before the first period.
    arma::vec v_h_init;

    /// n_factors; variance of the factor-innovation log-volatility innovations.
    arma::vec v_h_sigma;
};

/// Dynamic factor model with a normal prior on the loadings and on the factor
/// transition, and stochastic volatility in both error terms.
///
///     x_t = Lambda f_t + u_t,                  u_t ~ N(0, U_t),
///     f_t = sum_{j=1..p} A_j f_{t-j} + v_t,    v_t ~ N(0, V_t),
///
/// with U_t = diag(exp(h^u_t)) and V_t = diag(exp(h^v_t)), and each element of
/// both log-volatilities following a random walk of its own.
///
/// `Normal` names the coefficients, as it does in VarNormalStochvolInput against
/// VarTvpStochvolInput: the loadings and the factor transition are constant and
/// carry normal priors. Only the two error covariances move.
///
/// The two placements do different work and neither substitutes for the other.
/// Volatility in `u_t` reweights the series that identify the factors, so a
/// series that was noisy early and quiet later stops contributing on the same
/// terms throughout, and a single wild observation is absorbed where it happened
/// rather than dragged into the factor. Volatility in `v_t` is the common
/// component's own, and it is what stops the k idiosyncratic variances from
/// jointly absorbing a shock that every series felt at once -- the factor is
/// otherwise flattest exactly when it should move most.
///
/// Everything DfmNormalGammaInput says about the data and the two orderings
/// holds here unchanged: `train.y` is the only data, there is no ForecastData,
/// `spec.h` is the horizon, the free loadings run row by row, `a` is
/// vec([A_1 .. A_p]) with the blocks side by side, and the factors before the
/// sample are zero rather than drawn.
struct DfmNormalStochvolInput
{
    VarSpec spec;
    TrainData train;

    NormalPrior lambda_prior;  ///< Over the free loadings, in the row-major order above.
    NormalPrior a_prior;       ///< Unused when the factors have no dynamics.

    StochvolPrior u_sigma_prior;  ///< k idiosyncratic log-volatilities.
    StochvolPrior v_sigma_prior;  ///< n_factors factor-innovation log-volatilities.

    DfmNormalStochvolInitial initial;

    /// Whether the factors carry dynamics to draw. A transition order of zero is
    /// a static factor model with serially independent -- but still
    /// heteroskedastic -- factors, which this sampler estimates.
    bool use_a() const { return spec.n_factor_a() > 0; }

    void validate() const;
};

/// Where the DfmTvpGamma chain starts.
///
/// DfmNormalGammaInitial with the two coefficient blocks widened from a point to
/// a path, and each of them carrying the two extra pieces every random walk here
/// needs: the precision of its innovations and the state of the period before
/// the sample. The factor path is still not among these -- it is the first thing
/// every draw produces.
struct DfmTvpGammaInitial
{
    /// n_lambda x tt. Each column holds the free loadings of one period, row by
    /// row and within a row left to right -- the same ordering
    /// DfmNormalGammaInitial::lambda uses for its single vector, repeated once
    /// per period.
    arma::mat lambda;

    /// n_lambda x n_lambda; only the diagonal is read.
    arma::mat lambda_sigma_inv;

    /// n_lambda; the free loadings of the period before the sample.
    arma::vec lambda_init;

    /// n_factor_a x tt, each column vec([A_1 .. A_p]) of one period. Ignored
    /// when the factors have no dynamics.
    arma::mat a;

    /// n_factor_a x n_factor_a; only the diagonal is read.
    arma::mat a_sigma_inv;

    /// n_factor_a; the transition of the period before the sample.
    arma::vec a_init;

    /// k; the diagonal of the idiosyncratic precision, which does not move.
    arma::vec u_sigma_inv;

    /// n_factors; the diagonal of the factor innovation precision.
    arma::vec v_sigma_inv;
};

/// Dynamic factor model whose loadings and factor transition follow random
/// walks, with independent gamma priors on both error precisions.
///
///     x_t = Lambda_t f_t + u_t,                  u_t ~ N(0, U),  U diagonal,
///     f_t = sum_{j=1..p} A_{j,t} f_{t-j} + v_t,  v_t ~ N(0, V),  V diagonal,
///
/// with every free element of Lambda and every element of [A_1 .. A_p] a random
/// walk of its own. `Tvp` names the coefficients, exactly as it does in
/// VarTvpGammaInput against VarNormalGammaInput, and it names *both* coefficient
/// blocks: a model in which only the loadings drifted would be a different one,
/// and is not what this is.
///
/// The identifying block of Lambda does not drift, because it is not drawn: the
/// leading N x N block stays unit lower triangular in every period. Letting it
/// move would leave the rotation and scale of the factors free to wander over
/// the sample, and a loading path would then be a statement about the
/// normalisation as much as about the exposure it is read as.
///
/// Everything DfmNormalGammaInput says about the data holds here unchanged:
/// `train.y` holds the observed series and is the only data the model takes,
/// there is no ForecastData, `spec.h` is the horizon, the free loadings run row
/// by row, `a` is vec([A_1 .. A_p]) with the blocks side by side, and the
/// factors before the sample are zero rather than drawn.
///
/// Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
/// parameters: measuring changes in international business cycles. Federal
/// Reserve Bank of New York Staff Report No. 326.
struct DfmTvpGammaInput
{
    VarSpec spec;
    TrainData train;

    /// The state equation of the free loadings, in the row-major order above.
    RandomWalkPrior lambda_prior;

    /// The state equation of the factor transition. Unused when the factors
    /// have no dynamics.
    RandomWalkPrior a_prior;

    GammaPrior u_sigma_prior;  ///< k independent priors on the idiosyncratic precisions.
    GammaPrior v_sigma_prior;  ///< n_factors independent priors on the factor innovations.

    DfmTvpGammaInitial initial;

    /// Whether the factors carry dynamics to draw. A transition of order zero is
    /// a static factor model with serially independent factors, which this
    /// sampler estimates -- there is then no transition path either, and the
    /// model is a DFM with drifting loadings alone.
    bool use_a() const { return spec.n_factor_a() > 0; }

    void validate() const;
};

/// Where the DfmTvpStochvol chain starts.
///
/// DfmTvpGammaInitial with the two precisions replaced by a log-volatility path
/// each -- which is exactly what DfmNormalStochvolInitial does to
/// DfmNormalGammaInitial. Every one of the four blocks this model carries is a
/// path, and the factor path is still not among them: it is the first thing every
/// draw produces.
struct DfmTvpStochvolInitial
{
    /// n_lambda x tt, one column per period, each in the row-major order the
    /// sampler draws the free loadings in.
    arma::mat lambda;

    /// n_lambda x n_lambda; only the diagonal is read.
    arma::mat lambda_sigma_inv;

    /// n_lambda; the free loadings of the period before the sample.
    arma::vec lambda_init;

    /// n_factor_a x tt, each column vec([A_1 .. A_p]) of one period. Ignored
    /// when the factors have no dynamics.
    arma::mat a;

    /// n_factor_a x n_factor_a; only the diagonal is read.
    arma::mat a_sigma_inv;

    /// n_factor_a; the transition of the period before the sample.
    arma::vec a_init;

    /// tt x k log-volatilities of the idiosyncratic errors, one column per
    /// observed series. One period per *row*, unlike the coefficient paths
    /// above: that is the orientation the mixture routine works in, and the two
    /// are not interchangeable.
    arma::mat u_h;

    /// k; the idiosyncratic log-volatility before the first period.
    arma::vec u_h_init;

    /// k; variance of the idiosyncratic log-volatility innovations. A state
    /// rather than a prior -- the sampler redraws it every iteration -- even
    /// though the files keep it next to the prior it is drawn under.
    arma::vec u_h_sigma;

    /// tt x n_factors log-volatilities of the factor innovations.
    arma::mat v_h;

    /// n_factors; the factor-innovation log-volatility before the first period.
    arma::vec v_h_init;

    /// n_factors; variance of the factor-innovation log-volatility innovations.
    arma::vec v_h_sigma;
};

/// Dynamic factor model whose loadings and factor transition follow random walks
/// and whose two error terms carry stochastic volatility.
///
///     x_t = Lambda_t f_t + u_t,                  u_t ~ N(0, U_t),
///     f_t = sum_{j=1..p} A_{j,t} f_{t-j} + v_t,  v_t ~ N(0, V_t),
///
/// with U_t = diag(exp(h^u_t)) and V_t = diag(exp(h^v_t)), and every free element
/// of Lambda, every element of [A_1 .. A_p] and every element of both
/// log-volatilities a random walk of its own.
///
/// The widest model here: everything that can move does. It is DfmTvpGamma's
/// coefficients over DfmNormalStochvol's errors, and both halves earn their
/// place for the reasons those two set out -- a series' exposure to the common
/// component is not a constant of nature, and neither is the scale of the shocks
/// it is measured against. The two are also easy to confuse for one another in a
/// model that has only one of them: a series whose loading fell will look like a
/// series whose idiosyncratic variance rose, and a period of common turbulence
/// will look like a transition that changed. Carrying both is what lets the data
/// say which.
///
/// The identifying block does not drift, exactly as in DfmTvpGamma: Lambda's
/// leading N x N block stays unit lower triangular in every period, because it is
/// the normalisation rather than a parameter.
///
/// Everything DfmNormalGammaInput says about the data holds here unchanged:
/// `train.y` is the only data, there is no ForecastData, `spec.h` is the horizon,
/// the free loadings run row by row, `a` is vec([A_1 .. A_p]) with the blocks
/// side by side, and the factors before the sample are zero rather than drawn.
///
/// Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
/// parameters: measuring changes in international business cycles. Federal
/// Reserve Bank of New York Staff Report No. 326.
struct DfmTvpStochvolInput
{
    VarSpec spec;
    TrainData train;

    /// The state equation of the free loadings, in the row-major order above.
    RandomWalkPrior lambda_prior;

    /// The state equation of the factor transition. Unused when the factors
    /// have no dynamics.
    RandomWalkPrior a_prior;

    StochvolPrior u_sigma_prior;  ///< k idiosyncratic log-volatilities.
    StochvolPrior v_sigma_prior;  ///< n_factors factor-innovation log-volatilities.

    DfmTvpStochvolInitial initial;

    /// Whether the factors carry dynamics to draw. A transition of order zero is
    /// a static factor model with serially independent -- but still
    /// heteroskedastic -- factors and drifting loadings, which this sampler
    /// estimates.
    bool use_a() const { return spec.n_factor_a() > 0; }

    void validate() const;
};

/// Where the FavarNormalWishart chain starts.
///
/// The state path is not among these, for the reason DfmNormalGammaInitial
/// gives: it is the first thing every draw produces. Neither is the observed
/// half of the state, which is data and arrives in `train`.
struct FavarNormalWishartInitial
{
    /// The `n_favar_lambda` free loadings, in the order the sampler draws them
    /// -- row by row, and within a row left to right. See
    /// FavarNormalWishartInput.
    arma::vec lambda;

    /// vec([Phi_1 .. Phi_p]) of the state transition. Ignored when the state
    /// carries no dynamics.
    arma::vec a;

    /// k; the diagonal of the idiosyncratic precision. A vector rather than a
    /// matrix because R is diagonal by assumption -- that assumption is what
    /// makes this a factor model -- exactly as in DfmNormalGammaInitial.
    arma::vec u_sigma_inv;

    /// n_state x n_state precision of the state innovations. A matrix, not a
    /// diagonal, and that is the difference from the DFM beside it: see
    /// FavarNormalWishartInput.
    arma::mat v_sigma_inv;
};

/// Factor augmented VAR with a normal prior on the loadings and on the state
/// transition, independent gamma priors on the idiosyncratic precisions and a
/// Wishart prior on the precision of the state innovations.
///
///     x_t = Lambda_f f_t + Lambda_y y_t + e_t,   e_t ~ N(0, R),  R diagonal,
///     s_t = sum_{j=1..p} Phi_j s_{t-j} + v_t,    v_t ~ N(0, Q),
///
/// with s_t = (f_t', y_t')', for `k` observed series in the panel x_t,
/// `n_factors` unobserved factors f_t and `n_obs_factors` observed factors y_t,
/// after Bernanke, Boivin and Eliasz (2005).
///
/// What it is, against the dynamic factor models beside it. A DFM's state is
/// unobserved throughout; here half of it is data. The observed factors are the
/// variables whose dynamics the model is actually about -- a policy rate, output
/// -- and the panel is there to measure the common component they move with.
/// They sit in the state vector rather than beside it because the transition is
/// a VAR over both blocks jointly, and that coupling is the model: the factors
/// respond to the policy variable and the policy variable responds to the
/// factors.
///
/// Why the error precision is a Wishart here and a gamma in every DFM. R must
/// stay diagonal -- a factor model whose idiosyncratic errors may correlate has
/// nothing left for the factors to explain, which is why no `Dfm*` has the
/// option. Q is a different object: it is the innovation covariance of a VAR,
/// its observed block is an ordinary VAR covariance, and its cross block is the
/// correlation between the factor innovations and the shock to the observed
/// variables. Forcing that to zero would assert the policy shock is orthogonal
/// to every factor innovation, which is the one thing a FAVAR exists to
/// measure. So the third part of this model's name refers to Q, and R is
/// gamma-diagonal in every member of the family.
///
/// That choice is what fixes the identification, and the two cannot be picked
/// separately. A DFM pins the rotation of its factors with a unit lower
/// triangular loading block *and* a diagonal V. With Q free the second half is
/// gone, so the first half has to do the whole job: the leading `n_factors`
/// square block of Lambda_f is the *identity* here, and the observed columns of
/// those rows are zero. See VarSpec::n_favar_lambda(), which counts what that
/// leaves.
///
/// The data. `train.y` holds the panel, tt x k; `train.f_obs` holds the observed
/// factors, tt x n_obs_factors, one period per row. `z`, `x` and `w` are unused
/// -- a FAVAR has no regressors of the kind a VAR does. There is no ForecastData:
/// the forecast is the transition run forward, which needs no out-of-sample
/// matrix. `spec.h` is the horizon.
///
/// Two orderings have to be agreed with the host and are the only ones not
/// implied by a dimension:
///
///   - the free loadings, everywhere they appear as a vector -- `initial.lambda`
///     and both halves of `lambda_prior` -- run row by row from row `n_factors`
///     on, each row carrying its `n_state` elements with the factor loadings
///     before the observed ones. The rows before that are the identifying block
///     and carry none. That is the order the equation by equation draw consumes
///     them in. The *posterior* is stored differently, as vec of the whole
///     k x n_state matrix; see FavarNormalWishartDraws.
///   - `a` is vec([Phi_1 .. Phi_p]) with the blocks side by side, so the first
///     n_state^2 elements are Phi_1 read column-wise.
///
/// States before the sample are zero, both blocks of them, which is what makes
/// the first p transitions well defined -- the convention every DFM here
/// follows. It is a statement about the data: the observed factors are expected
/// in deviations, as a FAVAR's are anyway.
///
/// Bernanke, B. S., Boivin, J., & Eliasz, P. (2005). Measuring the effects of
/// monetary policy: a factor-augmented vector autoregressive (FAVAR) approach.
/// Quarterly Journal of Economics, 120(1), 387-422.
struct FavarNormalWishartInput
{
    VarSpec spec;
    TrainData train;

    NormalPrior lambda_prior;  ///< Over the free loadings, in the row-major order above.
    NormalPrior a_prior;       ///< Unused when the state has no dynamics.

    GammaPrior u_sigma_prior;    ///< k independent priors on the idiosyncratic precisions.
    WishartPrior v_sigma_prior;  ///< On the n_state square precision of the state innovations.

    FavarNormalWishartInitial initial;

    /// Whether the state carries dynamics to draw. An order of zero is a static
    /// model with serially independent states, which this sampler estimates --
    /// it is not an error, though it is a strange thing to ask a FAVAR for.
    bool use_a() const { return spec.n_favar_a() > 0; }

    /// Throws std::invalid_argument describing the first inconsistency it
    /// finds. Called by the sampler before the first draw, so a host that
    /// forgets to call it still gets the message rather than a crash a
    /// thousand iterations in.
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

/// Where the VecKlgs2010 chain starts. No inclusion indicators: the sampler
/// implements no variable selection.
struct VecKlgs2010Initial
{
    arma::vec a;
    arma::vec beta;
    arma::mat u_sigma_inv;
};

/// The complete argument of the VecKlgs2010 sampler -- VecNormalWishartInput
/// against the compact regressors.
///
/// `train.x` is read and `train.z` is not: this is the non-SUR reading of the
/// same model, so the regressors arrive one column per regressor rather than
/// kroneckered up with I_k. `forecast.z` is unaffected and stays in the level
/// VAR layout every VEC forecast expects, since the forecast is the level VAR's
/// either way.
struct VecKlgs2010Input
{
    VarSpec spec;
    TrainData train;
    ForecastData forecast; ///< `z` empty when no forecast was requested.

    NormalPrior a_prior;                ///< Unused when there are no regressors.
    ConstantCointSpacePrior beta_prior; ///< Unused when the rank is zero.
    WishartPrior u_sigma_prior;

    VecKlgs2010Initial initial;

    /// Coefficients in `a`, counted from the model dimensions rather than from
    /// the data.
    ///
    /// The other models read this off `z.n_cols`, which cannot be done here:
    /// the loadings have no column in `x` at all -- their regressor is
    /// beta' w_{t-1}, a function of the current draw rather than data -- so `x`
    /// is short by exactly the k*rank of them. Contemporaneous coefficients are
    /// left out because validate() refuses a structural model outright; a
    /// Wishart prior does not identify one.
    int n_a() const { return spec.n_alpha() + spec.n_non_structural_vec(); }

    /// Columns of the design matrix the sampler builds: the `rank` error
    /// correction columns in front of the compact regressors.
    int n_design() const { return spec.rank + spec.n_x_vec(); }

    /// Whether the model has coefficients to draw at all.
    bool use_a() const { return n_a() > 0; }
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
