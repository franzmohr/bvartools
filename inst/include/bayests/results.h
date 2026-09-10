// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_RESULTS_H
#define BAYESTS_RESULTS_H

#include "bayests/arma.h"

namespace bayests
{

/// Posterior draws of a VAR with a normal prior on the coefficients and a
/// Wishart prior on the error precision.
///
/// Draws run along the columns, which is the layout the samplers accumulate in
/// and the one that keeps a single draw contiguous. Hosts that want the
/// convention their ecosystem expects -- draws in rows, for both the HDF5
/// files and R -- transpose at the boundary.
struct VarNormalWishartDraws
{
    /// nparams x iterations. Empty when the model has no regressors.
    arma::mat a;

    /// nparams x iterations of zeros and ones. Empty unless variable
    /// selection was requested.
    arma::mat a_lambda;

    /// (k * k) x iterations; each column is a vectorised precision matrix.
    arma::mat u_sigma_inv;

    /// Length of the chain these draws came from.
    arma::uword iterations() const { return u_sigma_inv.n_cols; }

    bool has_a() const { return a.n_elem > 0; }
    bool has_lambda() const { return a_lambda.n_elem > 0; }
};

/// Posterior draws of a VAR with independent gamma priors on the error
/// precisions and, optionally, a constant covariance block.
///
/// Same convention as VarNormalWishartDraws: draws along the columns, and a
/// member left empty is how "the model did not have that part" is expressed.
struct VarNormalGammaDraws
{
    arma::mat a;           ///< nparams x iterations.
    arma::mat a_lambda;    ///< nparams x iterations of zeros and ones.

    /// (k * k) x iterations; each column a vectorised lower-triangular Psi.
    arma::mat psi;
    arma::mat psi_lambda;

    /// k x iterations. Only the diagonal is drawn, so only it is kept.
    arma::mat u_omega_inv;

    /// (k * k) x iterations; each column a vectorised precision matrix.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};

/// Posterior draws of a VAR with stochastic volatility.
///
/// The error precision moves with time, so unlike the constant-variance models
/// a draw is a whole path rather than a single matrix.
struct VarNormalStochvolDraws
{
    arma::mat a;        ///< nparams x iterations.
    arma::mat a_lambda;

    arma::mat psi;      ///< (k * k) x iterations.
    arma::mat psi_lambda;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};

/// Posterior draws of a VAR estimated at a conditional quantile.
///
/// The error variance moves with the latent scale of the asymmetric Laplace, so
/// like the stochastic volatility models a draw carries a whole path rather than
/// a single matrix. Unlike them there is no covariance block, so every stored
/// precision matrix is diagonal -- kept in the k x k form all the same, because
/// that is the shape every reader of this library expects and the one
/// iterations() is counted off.
///
/// The latent scales themselves are not kept. They are k * tt numbers per draw
/// of pure nuisance, and `u_omega_inv` together with `u_scale` recovers them.
struct VarNormalAldDraws
{
    arma::mat a;        ///< nparams x iterations.
    arma::mat a_lambda;

    /// k x iterations: the scale of the asymmetric Laplace, one per equation.
    arma::mat u_scale;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column. Diagonal throughout.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
};

/// Posterior draws of a VAR whose coefficients follow a random walk.
///
/// `a` holds a whole path per draw, which is what makes this model's output
/// shapes differ from the others': the coefficients are time-varying, and so,
/// when there is a covariance block, is Psi.
struct VarTvpGammaDraws
{
    arma::mat a;
    arma::mat a_sigma;
    arma::mat a_lambda;

    arma::mat psi;
    arma::mat psi_sigma;
    arma::mat psi_lambda;

    arma::mat u_omega_inv;
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};

/// Posterior draws of a VAR whose coefficients follow a random walk.
///
/// `a` holds a whole path per draw, which is what makes this model's output
/// shapes differ from the others': the coefficients are time-varying.
struct VarTvpWishartDraws
{
    arma::mat a;
    arma::mat a_sigma;
    arma::mat a_lambda;
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
};

/// Posterior draws of a VAR whose coefficients follow a random walk and whose
/// errors carry stochastic volatility.
///
/// Both the coefficients and the error precision are paths here, which is what
/// makes this the widest output of any of the models: `a` holds nparams * tt
/// per draw and `u_sigma_inv` holds k * k * tt.
struct VarTvpStochvolDraws
{
    arma::mat a;
    arma::mat a_sigma;
    arma::mat a_lambda;

    arma::mat psi;
    arma::mat psi_sigma;
    arma::mat psi_lambda;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};


/// Posterior draws of a VAR estimated at a conditional quantile whose
/// coefficients follow a random walk.
///
/// VarNormalAldDraws with the coefficients widened to a path and their state
/// innovation variance beside them, exactly as VarTvpStochvolDraws widens
/// VarNormalStochvolDraws.
struct VarTvpAldDraws
{
    arma::mat a;        ///< (nparams * tt) x iterations.
    arma::mat a_sigma;  ///< nparams x iterations.
    arma::mat a_lambda;

    /// k x iterations: the scale of the asymmetric Laplace, one per equation.
    arma::mat u_scale;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column. Diagonal throughout.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
};

/// Posterior draws of a VEC with a normal prior on the coefficients and a
/// Wishart prior on the error precision.
///
/// Draws run along the columns, which is the layout the samplers accumulate in
/// and the one that keeps a single draw contiguous. Hosts that want the
/// convention their ecosystem expects -- draws in rows, for both the HDF5
/// files and R -- transpose at the boundary.
struct VecNormalWishartDraws
{
    arma::mat a;          ///< nparams x iterations. Empty when the model has no regressors.
    arma::mat a_lambda;   ///< nparams x iterations of zeros and ones. Empty unless variable

    arma::mat beta;     ///< n_beta x iterations. Empty when the model has no cointegration.

    /// (k * k) x iterations; each column is a vectorised precision matrix.
    arma::mat u_sigma_inv;

    /// Length of the chain these draws came from.
    arma::uword iterations() const { return u_sigma_inv.n_cols; }

    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
    bool has_lambda() const { return a_lambda.n_elem > 0; }
};

/// Posterior draws of the non-SUR Koop, Leon-Gonzalez and Strachan (2010)
/// sampler.
///
/// The same posterior VecNormalWishartDraws carries, minus the inclusion
/// indicators: this sampler implements no variable selection. Laid out
/// identically, `a` holding vec(alpha) first, so a draw of one is a draw of the
/// other and the two can be compared coefficient by coefficient.
struct VecKlgs2010Draws
{
    arma::mat a;    ///< n_a x iterations. Empty when the model has no regressors.
    arma::mat beta; ///< n_beta x iterations. Empty without a cointegration relation.

    /// (k * k) x iterations; each column is a vectorised precision matrix.
    arma::mat u_sigma_inv;

    /// Length of the chain these draws came from.
    arma::uword iterations() const { return u_sigma_inv.n_cols; }

    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
};

/// Posterior draws of a VEC with independent gamma priors on the error
/// precisions and, optionally, a constant covariance block.
struct VecNormalGammaDraws
{
    arma::mat a;        ///< n_a x iterations.
    arma::mat a_lambda;

    arma::mat beta;     ///< n_beta x iterations. Empty without a cointegration relation.

    /// (k * k) x iterations; each column a vectorised lower-triangular Psi.
    arma::mat psi;
    arma::mat psi_lambda;

    /// k x iterations. Only the diagonal is drawn, so only it is kept.
    arma::mat u_omega_inv;

    /// (k * k) x iterations; each column a vectorised precision matrix.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
    bool has_lambda() const { return a_lambda.n_elem > 0; }
};

/// Posterior draws of a VEC with stochastic volatility.
///
/// The coefficients are constant and the error precision is not, so `a` and
/// `beta` hold one value per draw while `u_sigma_inv` holds a whole path.
struct VecNormalStochvolDraws
{
    arma::mat a;
    arma::mat a_lambda;

    arma::mat beta;

    arma::mat psi;      ///< (k * k) x iterations.
    arma::mat psi_lambda;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
    bool has_lambda() const { return a_lambda.n_elem > 0; }
};

/// Posterior draws of a VEC whose coefficients -- loadings and cointegration
/// vectors included -- follow a random walk, with a Wishart error precision.
struct VecTvpWishartDraws
{
    arma::mat a;        ///< (n_a * tt) x iterations.
    arma::mat a_sigma;
    arma::mat a_lambda;

    arma::mat beta;     ///< (n_beta * tt) x iterations.

    /// (k * k) x iterations; the precision does not move with time here.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
};

/// Posterior draws of a VEC whose coefficients -- loadings and cointegration
/// vectors included -- follow a random walk, with independent gamma priors on
/// the error precisions.
///
/// `u_sigma_inv` is (k * k) x iterations without a covariance block and
/// (k * k * tt) x iterations with one, since Psi is then a path and the
/// precision it implies moves with it. Same convention as VarTvpGamma.
struct VecTvpGammaDraws
{
    arma::mat a;
    arma::mat a_sigma;
    arma::mat a_lambda;

    arma::mat beta;

    arma::mat psi;
    arma::mat psi_sigma;
    arma::mat psi_lambda;

    arma::mat u_omega_inv; ///< k x iterations.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};

/// Posterior draws of a VEC whose coefficients -- loadings and cointegration
/// vectors included -- follow a random walk and whose errors carry stochastic
/// volatility.
///
/// Every member but the state variances is a path: `a` holds n_a * tt per draw,
/// `beta` holds n_beta * tt, and `u_sigma_inv` holds k * k * tt, periods stacked
/// within a column in each case.
struct VecTvpStochvolDraws
{
    arma::mat a;
    arma::mat a_sigma;
    arma::mat a_lambda;

    /// (n_beta * tt) x iterations. Empty when the model has no cointegration
    /// relation. `a` carries only the loadings on it, so this is the half
    /// without which Pi cannot be reconstructed.
    arma::mat beta;

    arma::mat psi;
    arma::mat psi_sigma;
    arma::mat psi_lambda;

    /// (k * tt) x iterations: the diagonal of the precision, period by period.
    arma::mat u_omega_inv;

    /// (k * k * tt) x iterations: one vectorised precision matrix per period,
    /// periods stacked within a column.
    arma::mat u_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_beta() const { return beta.n_elem > 0; }
    bool has_psi() const { return psi.n_elem > 0; }
};

/// Posterior draws of a dynamic factor model with independent gamma priors on
/// both error precisions.
///
/// Draws run along the columns, as everywhere else here. What is unusual about
/// this posterior is that the factors are part of it: they are unobserved
/// states, not parameters, so there is one whole path per draw and everything
/// downstream -- the pointwise log likelihood, the forecast -- reads it back
/// rather than re-filtering.
struct DfmNormalGammaDraws
{
    /// (k * n_factors) x iterations; each column is vec of the whole M x N
    /// loading matrix, the fixed ones and zeros of the identifying block
    /// included.
    ///
    /// Deliberately not the free elements alone, which is how the prior and the
    /// starting value are given. A caller wants Lambda, and reshaping k x
    /// n_factors gets it without having to know the identification rule; the
    /// free-element ordering is an implementation detail of the draw and stops
    /// at the edge of this struct.
    arma::mat lambda;

    /// (n_factors * tt) x iterations; each column is vec of the N x tt factor
    /// path, periods along the columns of that matrix.
    arma::mat factors;

    /// n_factor_a x iterations, vec([A_1 .. A_p]). Empty when the factors have
    /// no dynamics.
    arma::mat a;

    /// k x iterations. U is diagonal, so only the diagonal is drawn and only it
    /// is kept -- same convention as VarNormalGammaDraws::u_omega_inv.
    arma::mat u_sigma_inv;

    /// n_factors x iterations; the diagonal of the factor innovation precision.
    arma::mat v_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_factors() const { return factors.n_elem > 0; }
};

/// Posterior draws of a dynamic factor model with stochastic volatility in both
/// error terms.
///
/// DfmNormalGammaDraws with the two precisions widened from a point to a path.
/// The names are kept rather than renamed to the `u_omega_inv` the VAR and VEC
/// stochastic volatility models use for their diagonal: those models have both a
/// diagonal and a full precision to tell apart, and a dynamic factor model's
/// U and V are diagonal by assumption, so there is only one object and
/// DfmNormalGammaDraws already named it. A caller moving between the two DFMs
/// then finds the same fields, only wider.
struct DfmNormalStochvolDraws
{
    /// (k * n_factors) x iterations; each column is vec of the whole M x N
    /// loading matrix, the identifying block's fixed ones and zeros included.
    arma::mat lambda;

    /// (n_factors * tt) x iterations; each column is vec of the N x tt factor
    /// path, periods along the columns of that matrix.
    arma::mat factors;

    /// n_factor_a x iterations, vec([A_1 .. A_p]). Empty when the factors have
    /// no dynamics.
    arma::mat a;

    /// (k * tt) x iterations: the diagonal of the idiosyncratic precision,
    /// period by period, periods stacked within a column.
    arma::mat u_sigma_inv;

    /// (n_factors * tt) x iterations: the same for the factor innovations.
    arma::mat v_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_factors() const { return factors.n_elem > 0; }
};

/// Posterior draws of a dynamic factor model whose loadings and factor
/// transition follow random walks.
///
/// DfmNormalGammaDraws with the two coefficient blocks widened from a point to a
/// path, and a state variance added to each. The two precisions stay points:
/// this model's errors are homoskedastic, and it is the coefficients that move.
///
/// `lambda` and `a` are the two members whose height depends on which entry
/// point produced them, the same arrangement VarTvpGammaDraws has. The sampler
/// and the pointwise log likelihood carry the whole path -- every period under
/// its own coefficients -- while a forecast carries the last in-sample period
/// alone, since that is what it holds constant over the horizon.
struct DfmTvpGammaDraws
{
    /// (k * n_factors * tt) x iterations: one vectorised M x N loading matrix
    /// per period, periods stacked within a column, the identifying block's
    /// fixed ones and zeros included. Cut to k * n_factors for a forecast.
    ///
    /// Deliberately the whole matrix rather than the free elements the prior and
    /// the starting value are given as, for the reason DfmNormalGammaDraws gives:
    /// a caller wants Lambda_t, and reshaping k x n_factors gets it without
    /// having to know the identification rule.
    arma::mat lambda;

    /// n_lambda x iterations: the variance of the loading random walks, one per
    /// free element in the row-major order the prior uses.
    arma::mat lambda_sigma;

    /// (n_factors * tt) x iterations; each column is vec of the N x tt factor
    /// path, periods along the columns of that matrix.
    arma::mat factors;

    /// (n_factor_a * tt) x iterations, one vec([A_1 .. A_p]) per period. Empty
    /// when the factors have no dynamics; cut to n_factor_a for a forecast.
    arma::mat a;

    /// n_factor_a x iterations: the variance of the transition random walks.
    arma::mat a_sigma;

    /// k x iterations. U is diagonal, so only the diagonal is drawn and kept.
    arma::mat u_sigma_inv;

    /// n_factors x iterations; the diagonal of the factor innovation precision.
    arma::mat v_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_factors() const { return factors.n_elem > 0; }
};

/// Posterior draws of a dynamic factor model whose loadings, factor transition
/// and two error covariances all move with time.
///
/// DfmTvpGammaDraws with the two precisions widened from a point to a path, which
/// is exactly what DfmNormalStochvolDraws does to DfmNormalGammaDraws. Every
/// member but the two state variances is a path, which makes this the widest
/// output of any model here.
///
/// `lambda` and `a` are the members whose height depends on which entry point
/// produced them, the same arrangement VarTvpGammaDraws has, and `u_sigma_inv`
/// and `v_sigma_inv` are read the way DfmNormalStochvolDraws' are: whole for the
/// pointwise log likelihood, last period alone for a forecast.
struct DfmTvpStochvolDraws
{
    /// (k * n_factors * tt) x iterations: one vectorised M x N loading matrix
    /// per period, periods stacked within a column, the identifying block's
    /// fixed ones and zeros included. Cut to k * n_factors for a forecast.
    arma::mat lambda;

    /// n_lambda x iterations: the variance of the loading random walks, one per
    /// free element in the row-major order the prior uses.
    arma::mat lambda_sigma;

    /// (n_factors * tt) x iterations; each column is vec of the N x tt factor
    /// path, periods along the columns of that matrix.
    arma::mat factors;

    /// (n_factor_a * tt) x iterations, one vec([A_1 .. A_p]) per period. Empty
    /// when the factors have no dynamics; cut to n_factor_a for a forecast.
    arma::mat a;

    /// n_factor_a x iterations: the variance of the transition random walks.
    arma::mat a_sigma;

    /// (k * tt) x iterations: the diagonal of the idiosyncratic precision,
    /// period by period, periods stacked within a column.
    arma::mat u_sigma_inv;

    /// (n_factors * tt) x iterations: the same for the factor innovations.
    arma::mat v_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_factors() const { return factors.n_elem > 0; }
};

/// Posterior draws of a factor augmented VAR with a Wishart prior on the
/// precision of the state innovations.
///
/// DfmNormalGammaDraws with the state widened from the unobserved factors to
/// the unobserved and observed ones together, and the factor innovation
/// precision widened from a diagonal to a matrix. Draws run along the columns,
/// as everywhere else here.
///
/// Only the unobserved half of the state is a draw. The observed half is data
/// and is not copied into the posterior -- a caller that wants the whole state
/// path already holds it, and storing a second copy of the input once per draw
/// would be the largest thing in the file.
struct FavarNormalWishartDraws
{
    /// (k * n_state) x iterations; each column is vec of the whole
    /// k x n_state loading matrix [Lambda_f Lambda_y], the fixed ones and zeros
    /// of the identifying block included.
    ///
    /// Deliberately not the free elements alone, which is how the prior and the
    /// starting value are given, for the reason DfmNormalGammaDraws gives: a
    /// caller wants Lambda, and reshaping k x n_state gets it without having to
    /// know the identification rule.
    arma::mat lambda;

    /// (n_factors * tt) x iterations; each column is vec of the
    /// n_factors x tt path of the *unobserved* factors, periods along the
    /// columns of that matrix.
    arma::mat factors;

    /// n_favar_a x iterations, vec([Phi_1 .. Phi_p]) of the state transition.
    /// Empty when the state has no dynamics.
    arma::mat a;

    /// k x iterations. R is diagonal, so only the diagonal is drawn and only it
    /// is kept -- the convention DfmNormalGammaDraws set.
    arma::mat u_sigma_inv;

    /// (n_state * n_state) x iterations; each column is vec of the state
    /// innovation precision Q^-1. A whole matrix, unlike every DFM's, which is
    /// the point of this member of the family -- see FavarNormalWishartInput.
    arma::mat v_sigma_inv;

    arma::uword iterations() const { return u_sigma_inv.n_cols; }
    bool has_a() const { return a.n_elem > 0; }
    bool has_factors() const { return factors.n_elem > 0; }
};

/// Simulated forecast paths, (h * k) x draws: one column per posterior draw,
/// horizons stacked within a column in the same variable order as the sample.
///
/// A factor augmented VAR is the one model here whose forecast is wider than
/// `k`: it stacks (h * (k + n_obs_factors)), the panel of one horizon followed
/// by the observed factors of that horizon. The observed factors are what a
/// FAVAR is forecast for, and this is the only dataset that could carry them.
struct ForecastDraws
{
    arma::mat values;
};

} // namespace bayests

#endif // BAYESTS_RESULTS_H
