// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_SPEC_H
#define BAYESTS_SPEC_H

#include <string>

namespace bayests
{

/// Which variable selection scheme, if any, the sampler applies to the
/// coefficients in `a`.
enum class VarSelection
{
    none, ///< All regressors stay in the model.
    ssvs, ///< Stochastic search variable selection.
    bvs   ///< Bayesian variable selection.
};

/// Parses the spelling used on disk and at the R level. Throws
/// std::invalid_argument on anything else, rather than silently disabling
/// selection because a name was misspelled.
VarSelection var_selection_from_string(const std::string &name);

/// Inverse of var_selection_from_string.
const char *to_string(VarSelection selection);

/// The shape of the model, independent of the data and the priors.
struct VarSpec
{
    int k = 0; ///< Endogenous variables.
    int p = 0; ///< Lags of the endogenous variables.
    int m = 0; ///< Exogenous variables.
    int s = 0; ///< Lags of the exogenous variables.
    int n = 0; ///< Deterministic terms entering unrestricted.
    int h = 0; ///< Forecast horizon; 0 when no forecast was requested.

    /// Rank of the cointegration matrix. Zero for a VAR, and zero as well for
    /// a VEC estimated without a cointegration relation -- which is what makes
    /// the difference between a `beta` that is absent and one that is empty.
    int rank = 0;
    int k_beta = 0;

    /// Deterministic terms restricted to the cointegration space. Zero for a
    /// VAR. The unrestricted ones are counted by `n`, so the two never overlap
    /// and a VEC's deterministic terms are `n + n_restricted`.
    int n_restricted = 0;

    int iterations = 0; ///< Draws kept.
    int burnin = 0;     ///< Draws discarded before the first kept one.

    VarSelection varsel = VarSelection::none;

    /// Whether the error specification asks for a covariance block -- the
    /// "+covar" suffix the files carry. Which prefix it had, "gamma" or "sv",
    /// identifies the model and so is the reader's business, not the sampler's:
    /// by the time a sampler sees this, it already knows which one it is.
    bool covar = false;

    /// Structural model: the last k(k-1)/2 entries of `a` are contemporaneous
    /// coefficients rather than lag or deterministic terms, and are split off
    /// before a forecast path is simulated.
    bool structural = false;

    /// Total length of the chain.
    int draws() const { return iterations + burnin; }

    bool uses_varsel() const { return varsel != VarSelection::none; }

    /// Whether the model carries a Psi block. One variable has no off-diagonal
    /// covariance to model, so the flag on its own is not enough.
    bool uses_covar() const { return covar && k > 1; }

    /// Free elements of the lower triangle of Psi: k(k-1)/2, and zero unless
    /// the model has a covariance block.
    int n_psi() const { return uses_covar() ? k * (k - 1) / 2 : 0; }

    /// Contemporaneous coefficients carried at the end of `a`.
    int n_structural() const { return structural ? k * (k - 1) / 2 : 0; }

    /// Coefficients of a VAR that are not contemporaneous terms. Declared from
    /// the model's dimensions rather than counted off `z`, because a forecast
    /// runs from the posterior alone and has no `z` to count.
    ///
    /// A VAR regresses on p lags of the endogenous variables and on x_t together
    /// with s lags of it. A VEC does not -- see n_non_structural_vec().
    int n_non_structural() const { return k * (k * p + m * (s + 1) + n); }

    /// Whether the model carries a cointegration relation to draw. Zero rank
    /// leaves `beta` out of the posterior rather than storing an empty one.
    bool uses_coint() const { return rank > 0; }

    /// Free elements of `beta`: k_beta x rank, and zero unless the model has a
    /// cointegration relation.
    int n_beta() const { return k_beta * rank; }

    /// Loading coefficients carried at the front of `a`, and zero unless the
    /// model has a cointegration relation.
    int n_alpha() const { return k * rank; }

    /// Coefficients a VAR draw carries for one period, contemporaneous terms
    /// included. This is the width a stored coefficient path is cut into, so a
    /// host slicing one needs exactly this number and cannot recover it from the
    /// stored width alone.
    ///
    /// The loading term is included and is zero for every VAR, which is the only
    /// sense in which this is shared with a VEC: the *rest* of a VEC's count
    /// differs, so use nparams_per_period_vec() for one. Reading `p` and `s` as a
    /// VEC's own block counts is the mistake this pair of names exists to
    /// prevent.
    int nparams_per_period() const { return n_alpha() + n_non_structural() + n_structural(); }

    /// Lagged-difference blocks a VEC carries, Gamma_1..Gamma_{p-1}.
    ///
    /// `p` is the lag order of the *level* representation, the convention
    /// bvartools uses at the R level and the one the model files are written
    /// with, so differencing costs one block. A VEC of level order one has no
    /// Gamma at all.
    int n_gamma() const { return p > 1 ? p - 1 : 0; }

    /// Difference blocks of the unmodelled variables a VEC carries,
    /// Upsilon_0..Upsilon_{s-1}, the current difference dx_t included.
    ///
    /// `s` is a level order for the same reason `p` is, so there are s of these
    /// against the s + 1 level blocks x_t..x_{t-s} they imply.
    int n_upsilon() const { return m > 0 ? s : 0; }

    /// Regressors a VEC has besides its error correction term: the k(p-1)
    /// lagged differences of the endogenous variables, the m*s differences of
    /// the unmodelled ones and the n unrestricted deterministic terms.
    ///
    /// One column each, which is the compact per-period layout TrainData::x is
    /// written in and the width a k x n coefficient matrix has. The SUR layout
    /// spreads every one of them over k columns, so n_non_structural_vec() is k
    /// times this.
    int n_x_vec() const { return k * n_gamma() + m * n_upsilon() + n; }

    /// Coefficients of a VEC that are not loadings or contemporaneous terms:
    /// the Gamma blocks, the Upsilon blocks and the unrestricted deterministic
    /// terms. The restricted ones live in `beta`, not here.
    int n_non_structural_vec() const { return k * n_x_vec(); }

    /// Coefficients a VEC draw carries for one period, loadings and
    /// contemporaneous terms included -- the width of its `a`, and the number of
    /// columns its `z` has.
    int nparams_per_period_vec() const
    {
        return n_alpha() + n_non_structural_vec() + n_structural();
    }

    /// Throws std::invalid_argument if the counts cannot describe a model.
    /// Cheap, and the message is far better than the failure it prevents --
    /// a zero k turns the first reshape into a division by zero.
    void validate() const;
};

} // namespace bayests

#endif // BAYESTS_SPEC_H
