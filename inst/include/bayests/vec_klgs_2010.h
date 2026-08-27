// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_VEC_KLGS_2010_H
#define BAYESTS_VEC_KLGS_2010_H

#include "bayests/inputs.h"
#include "bayests/reporter.h"
#include "bayests/results.h"

namespace bayests
{

/// The Koop, Leon-Gonzalez and Strachan (2010) cointegration sampler written
/// against the compact regressors rather than the SUR system.
///
/// The model is VecNormalWishart's, exactly: a normal prior on the
/// non-cointegration coefficients, the cointegration space prior on beta and a
/// Wishart prior on the error precision, with the same three Gibbs blocks in the
/// same order and the same normalisation splitting each draw of Pi between alpha
/// and beta. What differs is only how the coefficient block is written down.
///
/// A VEC's k equations share their regressors, so its SUR design matrix is
/// z = kron(W_x, I_k), and the posterior precision the coefficient block needs
/// collapses accordingly:
///
///     z' kron(I_tt, Sigma^-1) z = kron(W_x' W_x, Sigma^-1),
///     z' kron(I_tt, Sigma^-1) y = vec(Sigma^-1 Y' W_x).
///
/// Both right-hand sides are what this sampler forms. The left-hand ones are
/// what VecNormalWishart forms, and the difference is not cosmetic: the SUR
/// spelling builds a (tt k) x (k n_x) regressor matrix and accumulates a Gram
/// product over it, which is O(tt k^3 n_x^2) against the O(tt n_x^2) of
/// W_x' W_x, and holds tt k^2 n_x doubles that the compact form never
/// materialises. For a six-variable model that is a factor of two hundred on the
/// hottest arithmetic in the chain.
///
/// The beta block is already compact in both -- it is linear in beta with
/// regressors kron(alpha, w_t'), which has no such structure to exploit -- so it
/// is the same code here as there, and so is Sigma's.
///
/// What is given up for it is variable selection. SSVS and BVS both act on
/// individual columns of the SUR design matrix, which is precisely the object
/// this sampler declines to build; validate() rejects either scheme and points
/// at VecNormalWishart, which draws the same posterior and carries both.
///
/// Values in, values out: no files, no console, no global state beyond the
/// Armadillo RNG. That is what lets the same object serve the command line and
/// an embedded caller such as an R package -- under RcppArmadillo the RNG is
/// R's own, so set.seed() reaches these draws without the sampler knowing.
///
/// Koop, G., Leon-Gonzalez, R., & Strachan, R. W. (2010). Efficient posterior
/// simulation for cointegrated models with priors on the cointegration space.
/// Econometric Reviews, 29(2), 224--242.
class VecKlgs2010Sampler
{
public:
    /// Runs the Gibbs sampler. Reports progress once per draw and honours an
    /// interrupt thrown from the reporter.
    ///
    /// Reads `input.train.x`, the compact regressors, and not `input.train.z`.
    ///
    /// Throws std::invalid_argument if `input` is inconsistent.
    VecKlgs2010Draws draw_coefficients(const VecKlgs2010Input &input, Reporter &reporter) const;

    /// Simulates one forecast path per posterior draw, in levels.
    ///
    /// The draws are rewritten in the level VAR parameterisation and the path is
    /// simulated by VarNormalWishartSampler, since a VEC and its level VAR are
    /// the same model and only one of them has a recursion worth writing twice.
    ///
    /// `input.forecast.z` is consequently expected in that level layout -- p + 1
    /// lags of the endogenous variables, s + 2 blocks of the unmodelled ones,
    /// then the unrestricted deterministic terms and the ones restricted to the
    /// cointegration space -- and not in the compact differenced layout
    /// `input.train.x` uses. The level history cannot be recovered from
    /// differenced regressors, so it has to be supplied rather than derived.
    ///
    /// Requires `input.forecast.z` and, when the model has coefficients,
    /// `draws.a`. Throws std::invalid_argument if either is missing, if spec.h
    /// is zero, or if the regressors do not match the converted coefficients.
    ForecastDraws forecast(const VecKlgs2010Input &input, const VecKlgs2010Draws &draws,
                           Reporter &reporter) const;

    /// Pointwise log likelihood, draws x periods -- one row per posterior
    /// draw, one column per observation, as expected by WAIC and PSIS-LOO.
    arma::mat log_likelihood(const VecKlgs2010Input &input, const VecKlgs2010Draws &draws) const;
};

} // namespace bayests

#endif // BAYESTS_VEC_KLGS_2010_H
