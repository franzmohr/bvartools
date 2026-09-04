// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_DATA_H
#define BAYESTS_DATA_H

#include "bayests/arma.h"

namespace bayests
{

/// The sample the model is estimated on.
struct TrainData
{
    /// Endogenous observations, one row per period and one column per
    /// variable. The samplers work with the stacked form vec(y'), so the
    /// column ordering is what fixes which coefficient in `a` belongs to which
    /// variable.
    ///
    /// A single row or a single column holding the already-stacked series is
    /// read the same way, which is how the HDF5 files store it.
    arma::mat y;

    /// Observations of the error correction term.
    arma::mat w;

    /// Regressors, (tt * k) rows by nparams columns, laid out to multiply the
    /// stacked response directly. Empty when the model has no regressors at
    /// all -- that is the switch the samplers test, not a separate flag.
    arma::mat z;

    /// The same regressors in the compact layout: one row per period and one
    /// column per regressor, so a coefficient matrix is k x n_x rather than a
    /// vector of length k * n_x.
    ///
    /// `z` is this kroneckered up with I_k, which is what lets the SUR samplers
    /// treat a multivariate model as one long regression -- at the price of a
    /// (tt k) x (k n_x) matrix and a Gram product k^3 times the size of this
    /// one's. A sampler reads whichever of the two it is written against and
    /// leaves the other empty; VecKlgs2010 is the one that reads this.
    arma::mat x;

    /// The observed factors of a factor augmented VAR, tt x n_obs_factors, one
    /// period per row. Empty for every other model.
    ///
    /// A member of its own rather than a reuse of `x`, which is the compact
    /// regressor layout and means something else. These are not regressors: they
    /// are the observed half of the state vector, they appear on the left of the
    /// transition as well as the right, and the model's own dynamics run over
    /// them. Handing them over as `x` would read as "a FAVAR is a regression
    /// with these on the right", which is the one thing it is not.
    arma::mat f_obs;

    /// Number of periods implied by `y` for a model with `k` variables.
    arma::uword periods(int k) const { return y.n_elem / static_cast<arma::uword>(k); }

    /// Number of coefficients in `a`.
    arma::uword nparams() const { return z.n_cols; }
};

/// The regressors the forecast is produced from.
struct ForecastData
{
    /// (h * k) rows by nparams columns. Lagged endogenous blocks are
    /// overwritten as the path is simulated, so only the deterministic and
    /// exogenous columns have to be filled in by the caller.
    arma::mat z;
};

} // namespace bayests

#endif // BAYESTS_DATA_H
