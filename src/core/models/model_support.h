// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_MODEL_SUPPORT_H
#define BAYESTS_CORE_MODELS_MODEL_SUPPORT_H

#include "bayests/data.h"
#include "bayests/priors.h"
#include "bayests/reporter.h"
#include "bayests/spec.h"
// Both unit lower triangular matrices a model carries are unpacked there, and
// the two are packed in different orders -- read that file before reaching for
// either. Included here so that a model needs one header, not two.
#include "core/algorithms/triangular_packing.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

namespace bayests::core
{

/// Tell the host, before the first draw, when BVS is about to select against a
/// prior it cannot select against.
///
/// `bayests check` reports the same thing from the same two functions, and a
/// run that is never checked would otherwise reach the end with inclusion
/// probabilities that describe the prior and nothing saying so. The Reporter is
/// how anything in here reaches a console: this file cannot print, and an R
/// package is not allowed to.
///
/// Called once per block, at the point the sweep's state is built, so a model
/// with a covariance block says it twice at most and a model selecting nothing
/// says nothing. Silent for every scheme but `bvs`, and for a prior tight
/// enough -- see bayests::flat_selection_prior(), which also says why only the
/// constant-coefficient models call this.
inline void report_flat_selection_prior(Reporter &reporter, const VarSelection scheme,
                                        const std::string &block, const VarSelPrior &prior,
                                        const arma::mat &v_inv)
{
    if (scheme != VarSelection::bvs)
    {
        return;
    }

    const FlatSelectionPrior report = flat_selection_prior(prior, v_inv);
    if (report.flat > 0)
    {
        reporter.message("warning: " + flat_selection_message(report, block));
    }
}

/// Refuses a NaN or an infinity anywhere in `values`, naming `what`.
///
/// BayesTS has no treatment of missing values, and a non-finite input does not
/// fail where it enters: it travels into a cross-product or a Cholesky factor
/// and surfaces as "inv_sympd(): matrix is singular" or "randg(): incorrect
/// distribution parameters", which names neither the dataset nor the value --
/// or, in a forecast or a score, it does not fail at all and comes back as a
/// NaN in the output of a run that reported success. Checked where the values
/// enter instead, so `bayests check` refuses the file the run would have.
///
/// The count rather than a position: an element's row and column read the
/// other way round in HDF5 dataspace terms than in Armadillo's, and a message
/// that names the wrong one is worse than one that names neither.
///
/// Read off the exponent bits, not std::isfinite(), for the reason all_finite()
/// in stochvol_mixture.h gives: a host compiling these sources with
/// -ffast-math is licensed to fold std::isfinite() to true, which would remove
/// exactly this check without a word.
inline void require_finite(const arma::mat &values, const std::string &what)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

    arma::uword bad = 0;
    for (arma::uword i = 0; i < values.n_elem; ++i)
    {
        const double value = values[i];
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        if ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL)
        {
            ++bad;
        }
    }
    if (bad > 0)
    {
        throw std::invalid_argument(
            what + " must be finite, but " + std::to_string(bad) + " of its " +
            std::to_string(values.n_elem) + " values " + (bad == 1 ? "is" : "are") +
            " NaN or infinite; missing values are not supported");
    }
}

/// The observations every model is given, checked with require_finite(): the
/// training sample and the realised values a forecast is scored against. Empty
/// members pass, so a model that reads only some of them needs no list of which.
///
/// Called from every input's validate(), and again where the score reads the
/// realised values -- scoring runs as a stage of its own, from draws a previous
/// run wrote, without validate() in between. The overload below adds the
/// forecast regressors, which the five factor models do not have: they
/// forecast from the factor transition alone.
inline void require_finite_observations(const TrainData &train, const TestData &test)
{
    require_finite(train.y, "the training observations /data/train/y");
    require_finite(train.z, "the regressors /data/train/z");
    require_finite(train.x, "the regressors /data/train/x");
    require_finite(train.w, "the cointegration regressors /data/train/w");
    require_finite(train.f_obs, "the observed factors /data/train/f_obs");
    require_finite(test.y, "the realised observations /data/test/y");
}

inline void require_finite_observations(const TrainData &train, const ForecastData &forecast,
                                        const TestData &test)
{
    require_finite_observations(train, test);
    require_finite(forecast.x, "the forecast regressors /data/forecast/x");
}

/// Rejects a forecast that was given no regressors by a model whose dimensions
/// say it has coefficients to apply to them.
///
/// `/data/forecast/x` is read with read_mat_if_present(), so a file that omits it
/// leaves this empty rather than failing. Empty then reads as "this model has no
/// regressors": use_a comes out false, the signal term drops out of the
/// recursion, and every horizon is drawn from the error distribution alone. That
/// path is written to /posterior/forecast/forecasts and reported as success, which
/// is worse than producing nothing -- nothing downstream can tell it from a model
/// whose coefficients are genuinely all zero.
///
/// Counted from the spec rather than from the posterior, so the message can say
/// what was expected even when the draws are missing as well.
inline void require_forecast_regressors(const VarSpec &spec, const arma::mat &x)
{
    if (spec.nparams_per_period() > 0 && x.n_elem == 0)
    {
        throw std::invalid_argument(
            "the model has " + std::to_string(spec.nparams_per_period()) +
            " coefficients but no forecast regressors were supplied; a forecast without them "
            "would be drawn from the error distribution alone");
    }
}

/// Rejects forecast regressors that do not have exactly one row per horizon.
///
/// Row i is read for horizon i and update_forecast_lags() writes the simulated
/// lags into it, so a matrix short of `h` rows is read and written past its end
/// -- an Armadillo exception in a checked build, and memory corruption in a
/// host that defines ARMA_NO_DEBUG. A longer one would run, on rows nothing
/// says the file meant to be ignored. `bayests check` refuses the same files.
///
/// An empty `x` passes: require_forecast_regressors() decides whether a model
/// with no regressors may forecast without them.
///
/// Finiteness is checked here as well, since this is where every VAR's and
/// VEC's forecast takes `x` and a forecast is a stage of its own that does not
/// go through validate(). The whole matrix, including the lag cells past the
/// first row that the forecast overwrites: a file is either free of missing
/// values or it is not, and a rule that depended on which cells a model reads
/// would be one more thing to get wrong.
inline void require_forecast_horizons(const arma::mat &x, const int h)
{
    require_finite(x, "the forecast regressors /data/forecast/x");
    if (x.n_elem > 0 && static_cast<arma::uword>(h) != x.n_rows)
    {
        throw std::invalid_argument("forecast regressors must have " + std::to_string(h) +
                                    " rows, one per horizon, got " +
                                    std::to_string(x.n_rows));
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

/// How far apart consecutive periods' precisions sit in a column of posterior
/// draws of u_sigma_inv: 0 where the column holds one k x k matrix for every
/// period, k^2 where it holds one per period, stacked.
///
/// The log likelihood of a model whose precision moves has to score each period
/// under that period's own. Reading the height rather than trusting a flag lets
/// one loop serve a model either way -- VarTvpGamma's precision moves exactly
/// when it has a covariance block -- and turns a draw of the wrong size into a
/// message rather than a reshape that silently keeps the first block.
inline arma::uword precision_stride(const arma::mat &u_sigma_inv, const int k, const int tt)
{
    const arma::uword kk = static_cast<arma::uword>(k) * k;
    if (u_sigma_inv.n_rows == kk)
    {
        return 0;
    }
    if (u_sigma_inv.n_rows == kk * tt)
    {
        return kk;
    }
    throw std::invalid_argument(
        "posterior draws of u_sigma_inv must have " + std::to_string(kk) +
        " rows, one matrix per draw, or " + std::to_string(kk * tt) + ", one per period, got " +
        std::to_string(u_sigma_inv.n_rows));
}

/// -log|Sigma| / 2 for the error covariance Sigma, read off its precision as
/// log|Sigma^-1| / 2 -- the determinant term of a Gaussian log likelihood.
///
/// Through a Cholesky in logs. The spelling this replaced inverted the precision
/// and took the log of the determinant of the result, which is a product of k
/// variances: it underflows to zero, and the log to minus infinity, once k and
/// the scale of the data are modest together -- twenty series with variances
/// near 1e-16 are enough. Symmetrised on the way in for the reason
/// draw_normal_precision() gives.
inline double half_log_det_precision(const arma::mat &precision)
{
    double value = 0.0;
    if (!arma::log_det_sympd(value, arma::mat(arma::symmatu(precision))))
    {
        throw std::runtime_error("a drawn error precision is not symmetric positive definite, so "
                                 "the log likelihood has no determinant term for it");
    }
    return value / 2;
}

/// tt identity matrices of order k, stacked row-wise: rows j k .. (j + 1) k - 1
/// are block j. The layout a time-varying Psi is held in, and the one the
/// per-period error precisions of the time-varying models share.
///
/// It replaces the (k tt) square block diagonal both used to be spelled as,
/// whose off-diagonal blocks were never read: 72 MB of zeros at k = 6,
/// tt = 500, and a dense product of order (k tt)^3 wherever Psi' Omega Psi was
/// formed from it rather than block by block.
inline arma::mat stacked_identity(const int k, const int tt)
{
    arma::mat stack(static_cast<arma::uword>(k) * tt, k, arma::fill::zeros);
    for (int j = 0; j < tt; j++)
    {
        stack.rows(j * k, (j + 1) * k - 1).diag().ones();
    }
    return stack;
}

/// The same fill, once per period, into the stack a time-varying Psi is held
/// as: rows j k .. (j + 1) k - 1 of `Psi` are the contemporaneous matrix of
/// period j, taken from column j of `psi`. See stacked_identity().
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
            Psi.submat(j * k + i, 0, j * k + i, i - 1) =
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
/// Writing the path in chronological order instead -- one write over
/// fcst[0 .. i*k-1], which is what every forecast here used to do -- reverses the
/// lags, putting A_1 on the oldest forecast rather than the newest. Only p <= 1
/// is insensitive to it, a single block having no order to get wrong, which is
/// why this survived: it needs p >= 2 and h >= 3 before the two spellings differ.
///
/// `x` is the compact layout, one period per row, so a lag block is k adjacent
/// entries of row i and the write is a copy. The SUR spelling this replaced put
/// the same k numbers through a kron with I_k and spread them over a k by k^2
/// submatrix, k^2 - k of whose entries were the zeros off that identity's
/// diagonal.
inline void update_forecast_lags(arma::mat &x, const arma::mat &fcst, const arma::uword draw,
                                 const int i, const int k, const int p)
{
    const int filled = i < p ? i : p;
    for (int j = 1; j <= filled; j++)
    {
        x.submat(i, (j - 1) * k, i, j * k - 1) =
            arma::trans(fcst.submat((i - j) * k, draw, (i - j + 1) * k - 1, draw));
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

/// The prior variance of a random walk's state before the sample, which is what
/// integrates that state out of the prior of the first period.
///
/// A random walk path is drawn by the simulation smoother, whose last two
/// arguments are the prior mean and covariance of the state the first observation
/// loads on. With a_1 = a_0 + v_1, v_1 ~ N(0, Sigma), and a_0 ~ N(mu_0, V_0), that
/// prior is N(mu_0, V_0 + Sigma) once a_0 is integrated out. Every time-varying
/// block passes exactly that, then draws a_0 given the path, and only then the
/// variance Sigma.
///
/// The blocks used to hand the smoother the a_0 of the previous draw and Sigma
/// itself, and to draw Sigma before a_0. Every one of those steps is a valid
/// conditional, but together they tie a_1 and a_0 to each other with variance
/// Sigma, and a random walk that is meant to move slowly has a small Sigma. With
/// a prior rate of 1e-12 on it the chain could not move the level of a path at
/// all: two chains started from different paths returned their own starting
/// values as the posterior, to the digits printed, and a time-varying VEC's
/// loadings stayed at whatever the host had initialised them to while its
/// cointegration vectors moved, which is how solved global models came out
/// explosive. Integrated out, a_0 leaves the level of the path to the data from
/// the first draw.
///
/// The order is part of the fix rather than a detail of it. a_0 was left out of
/// the path's draw, so it has to be drawn from its conditional on that path
/// before anything conditions on it again -- Sigma's innovations include
/// a_1 - a_0. Drawing Sigma first would make the sampler a partially collapsed
/// Gibbs sampler in an order that does not preserve the posterior (van Dyk and
/// Park, 2008).
///
/// Inverts the prior precision, which validate() therefore requires to be
/// positive definite for every time-varying block. It does not change over the
/// chain, so the callers take it once before the loop.
inline arma::mat initial_state_variance(const NormalPrior &prior)
{
    return arma::inv_sympd(prior.v_inv);
}

/// The two draws that follow a random walk coefficient path: the state of the
/// period before the sample, and the variance of its innovations. Both are
/// written in place, in that order.
///
/// `path` is n x tt, one period per column, as the simulation smoother returns
/// it, drawn with the state before the sample integrated out -- see
/// initial_state_variance(), which also says why `init` has to be drawn first.
/// `init` is normal, its only data being the first period of the path. The
/// innovations are then differences of the path against its own lag with the
/// period before the sample taken from the `init` just drawn, so their sum of
/// squares is what the inverse gamma posterior of the variance adds to the prior
/// rate.
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
/// time-varying models predate both and carry a copy each inline, in the same
/// order.
inline void draw_random_walk_state(arma::vec &sigma, arma::vec &init, const arma::mat &path,
                                   const arma::vec &post_shape, const arma::vec &prior_rate,
                                   const NormalPrior &init_prior)
{
    const arma::uword n = path.n_rows;
    const arma::uword tt = path.n_cols;

    const arma::mat init_precision = arma::diagmat(1.0 / sigma);
    init = draw_normal_precision(init_prior.v_inv + init_precision,
                                 init_prior.v_inv * init_prior.mu + init_precision * path.col(0));

    arma::mat differences(n, tt);
    differences.col(0) = path.col(0) - init;
    differences.cols(1, tt - 1) = path.cols(1, tt - 1) - path.cols(0, tt - 2);

    const arma::vec sse = arma::sum(arma::square(differences), 1);
    for (arma::uword i = 0; i < n; i++)
    {
        sigma(i) = 1.0 / arma::randg<double>(arma::distr_param(
                             post_shape(i), 1.0 / (prior_rate(i) + sse(i) * 0.5)));
    }
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
/// Called by `VarNormalStochvol` and the two stochastic volatility DFMs.
/// `VarTvpStochvol`, `VecNormalStochvol` and `VecTvpStochvol` still carry a copy
/// each, and `stochvol_mixture.h` says at length what came of the last pair of
/// copies in this library; they are left alone here only because rewriting a
/// sampler's draw sequence and rewriting this are separate changes.
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
