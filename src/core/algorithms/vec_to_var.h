// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_VEC_TO_VAR_H
#define BAYESTS_CORE_ALGORITHMS_VEC_TO_VAR_H

#include "bayests/results.h"
#include "bayests/spec.h"

namespace bayests::core
{

/// The VEC-to-VAR transformation: posterior draws of a vector error correction
/// model rewritten as the level VAR they imply.
///
/// A VEC and its level VAR are the same model in two parameterisations, so a
/// draw from one is a draw from the other -- there is nothing stochastic here,
/// only a change of basis applied draw by draw. What it buys is that everything
/// written for a VAR posterior (impulse responses, forecast error variance
/// decompositions, a forecast path in levels) applies to a VEC without being
/// written twice.
///
/// A VEC of level lag order `p` and rank `r` carries p - 1 lagged differences,
///
///     dy_t = alpha beta' w_{t-1} + sum_{i=1..p-1} Gamma_i dy_{t-i}
///            + sum_{j=0..s-1} Upsilon_j dx_{t-j} + C d_t + u_t,
///
/// with w_{t-1} = (y_{t-1}, x_{t-1}, d^r_{t-1}) and Pi = alpha beta' split
/// column-wise into Pi_y, Pi_x and Pi_d. The level VAR is then
///
///     y_t = sum_{i=1..p} A_i y_{t-i} + sum_{j=0..s} B_j x_{t-j}
///           + C d_t + Pi_d d^r_{t-1} + u_t,
///
///     A_1 = I + Pi_y + Gamma_1,  A_i = Gamma_i - Gamma_{i-1},  A_p = -Gamma_{p-1},
///     B_0 = Upsilon_0,  B_1 = Pi_x - Upsilon_0 + Upsilon_1,
///     B_j = Upsilon_j - Upsilon_{j-1},  B_s = -Upsilon_{s-1}.
///
/// Both recursions are the telescoping sum that substituting dy_t = y_t - y_{t-1}
/// into the VEC produces; see Lütkepohl (2006), New Introduction to Multiple
/// Time Series Analysis, section 6.3. Writing them out is equivalent to the
/// matrix product against the transformation matrices W and J that bvartools'
/// vec_to_var_transformation_matrix() builds, and costs O(k^2 p) rather than the
/// O(k^3 p^2) the product does.
///
/// How a VEC's VarSpec is read. `p` and `s` are the *level* orders, the same
/// convention bvartools uses at the R level, so they count one block more than
/// the VEC's own regressors do:
///
///   - `p`             level lags, hence Gamma_1..Gamma_{p-1}
///   - `m`, `s`        unmodelled variables and their level lags, hence
///                     Upsilon_0..Upsilon_{s-1}
///   - `n`             deterministic terms entering unrestricted, C
///   - `n_restricted`  deterministic terms restricted to the cointegration space
///   - `rank`          columns of alpha and beta
///
/// and `a` holds vec(alpha) first, then vec([Gamma_1 .. Gamma_{p-1}, Upsilon_0
/// .. Upsilon_{s-1}, C]), then the contemporaneous coefficients of a structural
/// model. Note that VarSpec::n_non_structural() and nparams_per_period() are the
/// VAR formulas and so do *not* describe a VEC's `a`: they would count k^2 p
/// where a VEC has k^2 (p - 1). This file counts the blocks itself rather than
/// asking them.

/// The dimensions of the level VAR implied by a VEC of shape `spec`.
///
/// `p` and `s` carry over unchanged, being level orders already; the restricted
/// deterministic terms fold into the unrestricted count, and no cointegration
/// relation is left to speak of. Variable selection is cleared with it: a level
/// coefficient is a sum of several VEC coefficients, so no single inclusion
/// indicator describes it -- see vec_to_var_coefficients().
///
/// The one adjustment is a floor of one lag. A VEC with p = 0 or p = 1 has no
/// Gamma at all, and its level form is the single block A_1 = I + Pi_y, so the
/// level order cannot be zero even when the differenced one is.
VarSpec vec_to_var_spec(const VarSpec &spec);

/// Rewrites `draws` in the level VAR representation, one draw at a time.
///
/// The result is laid out exactly as a VarNormalWishart posterior: `a` holds
/// vec([A_1 .. A_p, B_0 .. B_s, C, Pi_d]) followed by the
/// contemporaneous coefficients, which pass through untouched, and matches the
/// shape vec_to_var_spec(spec).nparams_per_period() reports. The deterministic
/// block keeps the unrestricted terms in the order the VEC had them and appends
/// the restricted ones, so a caller assembling the regressors for the level VAR
/// puts d_t first and the cointegration space's d^r_{t-1} last.
///
/// `u_sigma_inv` is shared, not transformed: the two parameterisations have the
/// same errors. `a_lambda` is left empty even when the VEC selected over its
/// coefficients, because inclusion does not survive the transformation -- A_i is
/// Gamma_i - Gamma_{i-1}, and one indicator cannot say which half was in.
///
/// Throws std::invalid_argument if `spec` cannot describe a model or if the
/// draws do not have the shape it describes. `spec.iterations` is not among the
/// things it insists on: the chain length comes from the posterior, so a
/// forecast-only input that never filled it in still converts.
VarNormalWishartDraws vec_to_var_coefficients(const VarSpec &spec,
                                              const VecNormalWishartDraws &draws);

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_VEC_TO_VAR_H
