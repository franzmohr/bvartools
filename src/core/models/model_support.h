// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_MODEL_SUPPORT_H
#define BAYESTS_CORE_MODELS_MODEL_SUPPORT_H

#include "bayests/data.h"
#include "bayests/priors.h"
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

/// The two draws that follow a random walk coefficient path: the variance of its
/// innovations, and the state of the period before the sample. Both are written
/// in place.
///
/// `path` is n x tt, one period per column, as the simulation smoother returns
/// it. The innovations are differences of the path against its own lag with the
/// period before the sample taken from `init`, so their sum of squares is what
/// the inverse gamma posterior of the variance adds to the prior rate; `init` is
/// then normal, its only data being the first period of the path.
///
/// `sigma` is the variance itself, not its inverse -- which is what the next draw
/// of the path is handed. `post_shape` is the prior shape plus tt/2, which does
/// not change over the chain and is formed once by the caller.
///
/// The same two draws draw_stochvol_state() makes below, over a path in the
/// other orientation: a log-volatility path is tt x k, one period per row,
/// because that is what the mixture routine works in. Sharing one function would
/// mean transposing one of them, which would reassociate the sums and move the
/// posteriors of four samplers to save a dozen lines. The VAR and VEC
/// time-varying models predate both and carry a copy each inline.
inline void draw_random_walk_state(arma::vec &sigma, arma::vec &init, const arma::mat &path,
                                   const arma::vec &post_shape, const arma::vec &prior_rate,
                                   const NormalPrior &init_prior)
{
    const arma::uword n = path.n_rows;
    const arma::uword tt = path.n_cols;

    arma::mat differences(n, tt);
    differences.col(0) = path.col(0) - init;
    differences.cols(1, tt - 1) = path.cols(1, tt - 1) - path.cols(0, tt - 2);

    const arma::vec sse = arma::sum(arma::square(differences), 1);
    for (arma::uword i = 0; i < n; i++)
    {
        sigma(i) = 1.0 / arma::randg<double>(arma::distr_param(
                             post_shape(i), 1.0 / (prior_rate(i) + sse(i) * 0.5)));
    }

    const arma::mat init_precision = arma::diagmat(1.0 / sigma);
    init = draw_normal_precision(init_prior.v_inv + init_precision,
                                 init_prior.v_inv * init_prior.mu + init_precision * path.col(0));
}

/// The two draws that follow the log-volatility path in a stochastic volatility
/// model: the variance of its random walk innovations, and the state of the
/// period before the sample. Both are written in place.
///
/// `h` is tt x K, one column per series, as `stochvol_ocsn_2007` returns it. The
/// random walk's innovations are h_t - h_{t-1} with h_0 taken from `h_init`, so
/// their sum of squares is what the inverse gamma posterior of the variance adds
/// to the prior rate. `h_init` is then normal, its only data being the first
/// period of the path -- every row of the random walk precision sums to zero but
/// the first.
///
/// `post_shape` is the prior shape plus tt/2, which does not change over the
/// chain and is formed once by the caller.
///
/// Factored out because a fifth model wanted it. The four that predate it carry
/// a copy each, and `stochvol_mixture.h` says at length what came of the last
/// pair of copies in this library; they are left alone here only because
/// rewriting a sampler's draw sequence and rewriting this are separate changes.
inline void draw_stochvol_state(arma::vec &h_sigma, arma::vec &h_init, const arma::mat &h,
                                const arma::vec &post_shape, const arma::vec &prior_rate,
                                const NormalPrior &h_init_prior)
{
    const arma::uword tt = h.n_rows;
    const arma::uword k = h.n_cols;

    arma::mat differences(tt, k);
    differences.row(0) = h.row(0) - arma::trans(h_init);
    differences.rows(1, tt - 1) = h.rows(1, tt - 1) - h.rows(0, tt - 2);

    const arma::vec sse = arma::trans(arma::sum(arma::square(differences)));
    for (arma::uword i = 0; i < k; i++)
    {
        h_sigma(i) = 1.0 / arma::randg<double>(arma::distr_param(
                               post_shape(i), 1.0 / (prior_rate(i) + sse(i) * 0.5)));
    }

    const arma::mat h_init_precision = arma::diagmat(1.0 / h_sigma);
    h_init = draw_normal_precision(h_init_prior.v_inv + h_init_precision,
                                   h_init_prior.v_inv * h_init_prior.mu +
                                       h_init_precision * arma::trans(h.row(0)));
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_MODEL_SUPPORT_H
