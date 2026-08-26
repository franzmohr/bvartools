// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_MODEL_SUPPORT_H
#define BAYESTS_CORE_MODELS_MODEL_SUPPORT_H

#include "bayests/data.h"
#include "bayests/spec.h"
// Both unit lower triangular matrices a model carries are unpacked there, and
// the two are packed in different orders -- read that file before reaching for
// either. Included here so that a model needs one header, not two.
#include "core/algorithms/triangular_packing.h"

#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Rejects a forecast that was given no regressors by a model whose dimensions
/// say it has coefficients to apply to them.
///
/// `/data/forecast/z` is read with read_mat_if_present(), so a file that omits it
/// leaves this empty rather than failing. Empty then reads as "this model has no
/// regressors": use_a comes out false, the signal term drops out of the
/// recursion, and every horizon is drawn from the error distribution alone. That
/// path is written to /posterior/forecast and reported as success, which is worse
/// than producing nothing -- nothing downstream can tell it from a model whose
/// coefficients are genuinely all zero.
///
/// Counted from the spec rather than from the posterior, so the message can say
/// what was expected even when the draws are missing as well.
inline void require_forecast_regressors(const VarSpec &spec, const arma::mat &z)
{
    if (spec.nparams_per_period() > 0 && z.n_elem == 0)
    {
        throw std::invalid_argument(
            "the model has " + std::to_string(spec.nparams_per_period()) +
            " coefficients but no forecast regressors were supplied; a forecast without them "
            "would be drawn from the error distribution alone");
    }
}

/// The response the samplers actually work with: the observations stacked
/// period by period, vec(y'). Storing `y` period-per-row and stacking here
/// keeps the caller's matrix in the orientation everyone else writes it in.
inline arma::vec stacked_response(const TrainData &train)
{
    return arma::vectorise(arma::trans(train.y));
}

/// One draw from the normal posterior written in precision form: returns a
/// sample from N(V^-1 b, V^-1), where `precision` is V and `rhs` is b, the
/// precision-weighted mean term.
///
/// Both the posterior mean and the draw need V factorised, and the obvious
/// spelling factorises it twice:
///
///     mu = arma::solve(V, b);                              // LU,   2n^3/3
///     x  = mu + arma::solve(arma::chol(V), arma::randn(n)); // chol,  n^3/3
///
/// `arma::solve()` on a general square matrix runs an LU with a reciprocal
/// condition estimate; it does not detect that V is symmetric positive
/// definite. Taking one Cholesky and reusing it for both the mean and the draw
/// costs n^3/3 in total -- a third of the work -- with the three triangular
/// solves being O(n^2) each.
///
/// V must be symmetric positive definite, which every posterior precision here
/// is by construction (prior precision plus a Gram matrix). `arma::chol()`
/// throws if it is not, which fails louder than the LU path did: an indefinite
/// V used to survive the mean solve and only fail on the draw.
inline arma::vec draw_normal_precision(const arma::mat &precision, const arma::vec &rhs)
{
    // Symmetrised on the way in. Every caller builds this as a prior precision
    // plus a Gram matrix, so it is symmetric in exact arithmetic -- but only in
    // exact arithmetic: z' D z is accumulated as two products, and the (i,j) and
    // (j,i) sums differ in their last bits. arma::chol() checks symmetry before
    // factorising and warns when the difference exceeds its tolerance, which a
    // flat prior makes easy to reach because there is nothing on the diagonal to
    // dominate it. Reflecting the upper triangle costs one pass and is exactly
    // right for a matrix that is symmetric up to rounding; the alternative is a
    // warning on every draw and, when the asymmetry grows, a failure.
    const arma::mat r = arma::chol(arma::symmatu(precision));

    // mean = precision^-1 rhs, by forward then back substitution.
    const arma::vec mean = arma::solve(arma::trimatu(r),
                                       arma::solve(arma::trimatl(r.t()), rhs));

    // Cov(r^-1 z) = r^-1 r^-T = (r.t() * r)^-1 = precision^-1.
    return mean + arma::solve(arma::trimatu(r), arma::randn<arma::vec>(precision.n_rows));
}

/// The same fill, once per period, into the block diagonal a time-varying Psi
/// is stored as: block j of `Psi` is the contemporaneous matrix of period j,
/// taken from column j of `psi`.
///
/// Not expressed in terms of fill_strict_lower_triangle(): the destination is a
/// submatrix of `Psi` rather than a matrix, and an Armadillo subview cannot be
/// sliced again.
inline void fill_psi_path(arma::mat &Psi, const arma::mat &psi, const int k)
{
    const int tt = static_cast<int>(psi.n_cols);
    for (int j = 0; j < tt; j++)
    {
        for (int i = 1; i < k; i++)
        {
            Psi.submat(j * k + i, j * k, j * k + i, j * k + i - 1) =
                arma::trans(psi.submat(i * (i - 1) / 2, j, (i + 1) * i / 2 - 1, j));
        }
    }
}

/// Splits the contemporaneous coefficients off the end of `a`, returning them
/// and shortening `a` to the coefficients that do have a column in `z`.
///
/// A structural model carries its k(k-1)/2 contemporaneous coefficients as the
/// last rows of the posterior, and they have no regressors: `z.n_cols` is short
/// by exactly that many. So `nparams` has to be the *posterior's* own count,
/// which the callers derive two different and equally correct ways -- off
/// `coefficients.a` where the coefficients are constant, off the spec where they
/// are a path and the posterior holds one period. Splitting on `z.n_cols`
/// instead cuts `a` in the wrong place: it takes the contemporaneous block out
/// of the lag coefficients and leaves a width that no longer matches `z`.
///
/// Returns an empty matrix for a model that is not structural, which is then the
/// flag the caller tests -- there is nothing to split and nothing to apply.
inline arma::mat split_structural_coefficients(const VarSpec &spec, arma::mat &a,
                                               const int nparams)
{
    if (!spec.structural)
    {
        return {};
    }

    const int n_structural = spec.n_structural();
    arma::mat a0 = a.rows(nparams - n_structural, nparams - 1);

    // A model that is nothing but its contemporaneous coefficients leaves `a`
    // alone: there is no row left to keep, and the caller's use_a is false.
    if (nparams > n_structural)
    {
        a = a.rows(0, nparams - n_structural - 1);
    }

    return a0;
}

/// A_0^{-1} for one draw, unpacked from the block split off above.
///
/// By column, unlike Psi -- see core/algorithms/triangular_packing.h, which is
/// where the two orders and the reason they differ are written down.
///
/// Called once per draw. The two samplers that had this written out inline did
/// it inside the horizon loop instead, rebuilding and re-inverting the same
/// matrix h times; nothing in it depends on the horizon.
inline arma::mat structural_inverse(const arma::mat &a0, const arma::uword draw,
                                    const arma::mat &diag_k)
{
    arma::mat a_0 = diag_k;
    fill_strict_lower_triangle_by_column(a_0, a0.col(draw));
    return arma::solve(a_0, diag_k);
}

/// Writes the simulated path into the lagged-endogenous columns of a forecast's
/// regressor matrix, for horizon `i` of draw `draw`.
///
/// The lag blocks run most recent first: column block j carries y_{t-j}, which is
/// how the training regressors are laid out -- verified against the recorded
/// fixtures, whose first observation has the immediately preceding period in
/// block one and the one before it in block two. At horizon i the block for lag j
/// is therefore the forecast made for horizon i - j, and the blocks past lag i
/// are still actual observations, which the caller supplied and this leaves alone.
///
/// Writing the path in chronological order instead -- one kron over
/// fcst[0 .. i*k-1], which is what every forecast here used to do -- reverses the
/// lags, putting A_1 on the oldest forecast rather than the newest. Only p <= 1
/// is insensitive to it, a single block having no order to get wrong, which is
/// why this survived: it needs p >= 2 and h >= 3 before the two spellings differ.
inline void update_forecast_lags(arma::mat &z, const arma::mat &fcst, const arma::uword draw,
                                 const int i, const int k, const int p, const arma::mat &diag_k)
{
    const int filled = i < p ? i : p;
    for (int j = 1; j <= filled; j++)
    {
        z.submat(i * k, (j - 1) * k * k, (i + 1) * k - 1, j * k * k - 1) =
            arma::kron(arma::trans(fcst.submat((i - j) * k, draw, (i - j + 1) * k - 1, draw)),
                       diag_k);
    }
}

/// The regressors of the psi block: equation i of period j is explained by the
/// errors of the equations above it, so row j(k-1)+i-1 carries -u(0..i-1, j) in
/// the columns belonging to row i of Psi.
///
/// `u` is the k x tt error matrix and `psi_z` is expected at its full size,
/// tt(k-1) x k(k-1)/2, with the cells outside those blocks left alone -- they
/// are structurally zero and stay zero for the life of the chain.
inline void build_psi_regressors(arma::mat &psi_z, const arma::mat &u)
{
    const int k = static_cast<int>(u.n_rows);
    const int tt = static_cast<int>(u.n_cols);
    for (int i = 1; i < k; i++)
    {
        for (int j = 0; j < tt; j++)
        {
            psi_z.submat(j * (k - 1) + i - 1, i * (i - 1) / 2,
                         j * (k - 1) + i - 1, (i + 1) * i / 2 - 1) =
                -arma::trans(u.submat(0, j, i - 1, j));
        }
    }
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_MODEL_SUPPORT_H
