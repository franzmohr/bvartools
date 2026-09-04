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

    /// Unobserved factors a dynamic factor model carries, the N of
    ///
    ///     x_t = Lambda f_t + u_t,   f_t = sum_j A_j f_{t-j} + v_t.
    ///
    /// Zero for every model that is not a factor model, exactly as `rank` is
    /// zero for every model that is not a VEC. A factor model reads the rest of
    /// this struct its own way and only two other fields carry a meaning for
    /// one: `k` is the number of *observed* series, the M above, and `p` is the
    /// lag order of the factor transition. `m`, `s`, `n`, `rank` and
    /// `n_restricted` are all zero.
    int n_factors = 0;

    /// Observed factors a factor augmented VAR carries, the Ky of
    ///
    ///     x_t = Lambda_f f_t + Lambda_y y_t + e_t,
    ///     s_t = sum_j Phi_j s_{t-j} + v_t,   s_t = (f_t', y_t')'.
    ///
    /// Zero for every model that is not a FAVAR, a dynamic factor model
    /// included -- which is exactly what a FAVAR with none of them would be.
    ///
    /// `y_t` is data and not a state to draw, and it sits in the state vector
    /// regardless: the transition is a VAR over the two blocks jointly, and that
    /// coupling is the whole of what separates this from the DFM beside it. It
    /// is why a FAVAR's state path draw conditions rather than draws, and why
    /// `n_state()` rather than `n_factors` is the width of everything the
    /// transition touches.
    int n_obs_factors = 0;

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

    /// Whether the model carries unobserved factors to draw.
    bool uses_factors() const { return n_factors > 0; }

    /// Free elements of a DFM's M x N loading matrix.
    ///
    /// Lambda is not free throughout: only the product Lambda f_t is identified,
    /// so its leading N x N block is fixed unit lower triangular -- ones on the
    /// diagonal, zeros above -- which pins both the rotation and the scale of the
    /// factors. That leaves min(i, N) free elements in row i, and summing them
    /// over the M rows gives N(2M - N - 1)/2.
    ///
    /// Zero unless the model has factors, and zero as well if `n_factors`
    /// exceeds `k`, which is not a model -- validate() rejects it, and this
    /// declines to return a negative count in the meantime.
    int n_lambda() const
    {
        return (uses_factors() && n_factors <= k) ? n_factors * (2 * k - n_factors - 1) / 2 : 0;
    }

    /// Coefficients of a DFM's factor transition, vec([A_1 .. A_p]). Zero when
    /// the factors carry no dynamics of their own.
    int n_factor_a() const { return n_factors * n_factors * p; }

    /// Whether the model carries observed factors beside the unobserved ones.
    bool uses_obs_factors() const { return n_obs_factors > 0; }

    /// Elements of the state vector a FAVAR's transition runs over: the
    /// unobserved factors and the observed ones together. Equal to `n_factors`
    /// for a DFM, which is the case of none of the second kind.
    int n_state() const { return n_factors + n_obs_factors; }

    /// Free elements of a FAVAR's k x n_state loading matrix: the whole width of
    /// every row from `n_factors` on, and nothing at all before that, so
    /// (k - n_factors) n_state.
    ///
    /// The leading block is fixed, as in a DFM, and it is fixed *harder*: the
    /// n_factors square block of the factor columns is the identity rather than
    /// a unit lower triangle, and the observed columns of those same rows are
    /// zero. The first n_factors series of the panel are therefore the factors
    /// plus idiosyncratic noise and nothing else.
    ///
    /// That is not a stylistic difference from n_lambda() and the two do not
    /// coincide at any dimension. A rotation F -> C F leaves the measurement
    /// alone if Lambda_f absorbs it, so something has to rule C out. A dynamic
    /// factor model rules it out with two restrictions at once -- a unit lower
    /// triangular block and a *diagonal* V, which together admit only C = I,
    /// the uniqueness of an LDL factorisation. A FAVAR cannot use the second
    /// half of that: its Q is unrestricted by design, that being the whole point
    /// of the model. A unit lower triangular block alone leaves every unit lower
    /// triangular C admissible, and the model is then not identified -- it runs,
    /// it produces plausible numbers, and its loadings wander along a ridge. The
    /// identity block rules out C on its own, which is why Bernanke, Boivin and
    /// Eliasz use it and why this does.
    ///
    /// The two restrictions cost the same. A unit triangle leaves
    /// n_factors(n_factors - 1)/2 loadings free and spends that many on making V
    /// diagonal; the identity spends them on the loadings and leaves Q free.
    ///
    /// Paired with n_lambda() rather than replacing it, the way
    /// n_non_structural_vec() is paired with n_non_structural(): reading a
    /// FAVAR's loadings off the DFM's count is the mistake the pair exists to
    /// prevent, and here it is a mistake at every dimension rather than only
    /// when there are observed factors.
    ///
    /// Zero unless the model has factors, and zero as well if `n_factors`
    /// exceeds `k`, on the same grounds as n_lambda().
    int n_favar_lambda() const
    {
        return (uses_factors() && n_factors <= k) ? (k - n_factors) * n_state() : 0;
    }

    /// Coefficients of a FAVAR's state transition, vec([Phi_1 .. Phi_p]). The
    /// same pairing with n_factor_a(), and the same reason: the transition is a
    /// VAR over the whole state, not over the factors alone.
    int n_favar_a() const { return n_state() * n_state() * p; }

    /// Throws std::invalid_argument if the counts cannot describe a model.
    /// Cheap, and the message is far better than the failure it prevents --
    /// a zero k turns the first reshape into a division by zero.
    void validate() const;
};

} // namespace bayests

#endif // BAYESTS_SPEC_H
