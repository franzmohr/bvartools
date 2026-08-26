// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_TRIANGULAR_PACKING_H
#define BAYESTS_CORE_ALGORITHMS_TRIANGULAR_PACKING_H

#include "bayests/arma.h"

namespace bayests::core
{

/// Unpacking the two unit lower triangular matrices a model can carry.
///
/// Both Psi and the structural matrix A_0 have ones on the diagonal and nothing
/// above it, so both are stored as their k(k-1)/2 free elements alone -- and the
/// two are stored in *different orders*, which is the reason this file exists
/// rather than one function serving both:
///
///   Psi   row by row:      (1,0) (2,0) (2,1) (3,0) (3,1) (3,2) ...
///   A_0   column by column: (1,0) (2,0) (3,0) (2,1) (3,1) (3,2) ...
///
/// They agree up to k = 3 and diverge from k = 4 on, which is exactly the shape
/// of bug this pair of names is meant to prevent: three variables is the size
/// most tests and fixtures use, so picking the wrong one is invisible until
/// someone estimates a four-variable model.
///
/// The orders are not a choice made here. They are the layouts the regressors
/// arrive in, and each was read off the matrix that produces them:
///
///   * Psi's regressors are built by build_psi_regressors(), which walks
///     equations outermost -- row i of Psi is explained by the errors of the
///     i equations above it, contiguously.
///   * A_0's regressors are the columns of kron(-y, I_k) that survive dropping
///     the diagonal and everything above it. That kron puts variable j's block
///     of k columns together, one per equation, so deleting the upper triangle
///     leaves the strict lower one ordered by variable -- by column.
///
/// This file lives in algorithms/ rather than beside the models that use it
/// because vec_to_var.cpp needs the second one too, and an algorithm reaching
/// up into models/ would invert the layering.

/// Writes `packed` into the strict lower triangle of `dest`, row by row: row i
/// holds i elements, starting at i(i-1)/2. This is Psi's order -- for A_0 see
/// fill_strict_lower_triangle_by_column().
///
/// `dest` must already hold the unit diagonal; only the strict lower triangle is
/// touched. Works whether `dest` is dense or sparse and whether the packed
/// source is a draw or a vector of inclusion indicators.
template <class Matrix, class Packed>
inline void fill_strict_lower_triangle(Matrix &dest, const Packed &packed)
{
    for (arma::uword i = 1; i < dest.n_rows; i++)
    {
        dest.submat(i, 0, i, i - 1) =
            arma::trans(packed.subvec(i * (i - 1) / 2, (i + 1) * i / 2 - 1));
    }
}

/// Writes `packed` into the strict lower triangle of `dest`, column by column:
/// column j holds the k - 1 - j elements below the diagonal. This is A_0's
/// order -- for Psi see fill_strict_lower_triangle().
///
/// Element-at-a-time rather than a submatrix assignment per column: the source
/// is contiguous but the destination is a column of a column-major matrix, and
/// an Armadillo subview of a subview cannot be assigned to.
///
/// `dest` must already hold the unit diagonal; only the strict lower triangle is
/// touched.
template <class Matrix, class Packed>
inline void fill_strict_lower_triangle_by_column(Matrix &dest, const Packed &packed)
{
    arma::uword at = 0;
    for (arma::uword j = 0; j + 1 < dest.n_cols; j++)
    {
        for (arma::uword i = j + 1; i < dest.n_rows; i++)
        {
            dest(i, j) = packed(at++);
        }
    }
}

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_TRIANGULAR_PACKING_H
