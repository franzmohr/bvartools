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

/// Simulated forecast paths, (h * k) x draws: one column per posterior draw,
/// horizons stacked within a column in the same variable order as the sample.
struct ForecastDraws
{
    arma::mat values;
};

} // namespace bayests

#endif // BAYESTS_RESULTS_H
