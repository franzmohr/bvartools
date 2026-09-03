// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_ALGORITHMS_STOCHVOL_MIXTURE_H
#define BAYESTS_CORE_ALGORITHMS_STOCHVOL_MIXTURE_H

#include "bayests/arma.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

/**
 * @file stochvol_mixture.h
 * @brief The stochastic volatility draw that both mixture approximations share.
 *
 * Kim, Shephard and Chib (1998) and Omori, Chib, Shephard and Nakajima (2007)
 * are the same algorithm run with a different table: seven normal components
 * against ten, approximating the same log chi-squared distribution. Everything
 * around the table -- the checks, the indicator draw, the banded draw of the
 * path -- is identical, so it lives here and each of the two translation units
 * next door is its constants plus a call.
 *
 * That is not a tidiness argument. The two were copies, and the copy is what
 * produced the one bug either of them has had: the component count of the
 * ten-component mixture was carried over into the seven-component one, which
 * left three tables padded with zeros, divided by a zero variance and returned a
 * singular draw on the first call. A fix applied to one copy and not the other
 * is the same failure waiting, and the log-domain indicator draw below had
 * already been through exactly that.
 */

namespace bayests::core
{

/// A normal mixture approximating the log chi-squared distribution of
/// \f$\log(u^2)\f$, which is what turns the non-linear measurement equation of a
/// stochastic volatility model into a conditionally linear one.
///
/// The three tables are `arma::rowvec` rather than `arma::rowvec::fixed<N>` on
/// purpose: a plain `rowvec` takes its length from the initialiser list, so a
/// table cannot end up longer than the values written into it and the count
/// cannot disagree with the tables. `-DARMA_NO_DEBUG`, which an embedded host
/// sets, removes the size check that would otherwise catch it.
///
/// The means are the published ones less 1.2704, the mean of the log
/// chi-squared distribution, so that the mixture approximates it centred.
struct NormalMixture
{
    arma::rowvec weight;
    arma::rowvec mean;
    arma::rowvec variance;

    arma::uword size() const
    {
        return weight.n_elem;
    }
};

namespace stochvol_detail
{

/// Everything the caller can get wrong about an argument is one message, so
/// that a failure names the argument rather than surfacing as an Armadillo
/// bounds error or a matrix of NaNs several draws later.
inline void require(const char *algorithm, bool ok, const std::string &what)
{
    if (!ok)
    {
        throw std::invalid_argument(std::string(algorithm) + ": " + what);
    }
}

inline std::string dims(const arma::mat &m)
{
    return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

/// True if every element is a finite number, read off the exponent bits.
///
/// Not `arma::is_finite()` and not `std::isfinite()`: a build with
/// `-ffast-math` gets `-ffinite-math-only` with it, which licenses the compiler
/// to fold either of them to `true`. Armadillo says as much, with a warning on
/// stderr per call that a sampler would emit once per draw. This project does
/// not set the flag, and the top-level CMakeLists says why, but an embedded host
/// compiles these sources under its own toolchain's rules. Looking at the bit
/// pattern instead is an integer operation that no floating point assumption can
/// remove, and an all-ones exponent is exactly an infinity or a NaN in IEEE-754
/// binary64.
inline bool all_finite(const arma::mat &m)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

    const arma::uword n = m.n_elem;
    for (arma::uword i = 0; i < n; i++)
    {
        std::uint64_t bits;
        const double value = m[i];
        std::memcpy(&bits, &value, sizeof(bits));
        if ((bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL)
        {
            return false;
        }
    }
    return true;
}

inline void require_length(const char *algorithm, const arma::vec &v, arma::uword n,
                           const char *what)
{
    require(algorithm, v.n_elem == n,
            std::string(what) + " must have " + std::to_string(n) + " elements, got " +
                std::to_string(v.n_elem));
}

/// Every element has to be finite and strictly positive: `sigma` is divided by
/// and `constant` is added under a logarithm, so a zero or a NaN here becomes an
/// infinity that only shows up as a broken draw further down.
inline void require_positive(const char *algorithm, const arma::vec &v, const char *what)
{
    require(algorithm, all_finite(v), std::string(what) + " contains NaN or infinite values");
    require(algorithm, v.min() > 0.0,
            std::string(what) + " must be strictly positive, but element " +
                std::to_string(v.index_min() + 1) + " is " + std::to_string(v.min()));
}

/// Draws the mixture indicator of every period for one variable.
///
/// The component probabilities are formed in logs and shifted by their row
/// maximum before they are exponentiated. The direct spelling -- weight times
/// `arma::normpdf`, then divide by the row sum -- underflows to a row of zeros
/// as soon as an observation sits far enough out in the tails of every
/// component, and the normalisation then turns the row into NaNs and the draw
/// into an out-of-range index. After the shift the largest entry of every row is
/// exactly one, so the row sum can never be zero.
///
/// @param y_star T-vector of the transformed observations \f$\log(u_t^2 + c)\f$.
/// @param h_col T-vector of the current log-volatility.
/// @param mixture the components to choose among.
/// @return T-vector of indices in `[0, mixture.size() - 1]`.
inline arma::uvec draw_mixture_states(const arma::vec &y_star, const arma::vec &h_col,
                                      const NormalMixture &mixture)
{
    const arma::uword tt = y_star.n_elem;
    const arma::uword n = mixture.size();

    // log p_j - 0.5 log sigma2_j - 0.5 (y*_t - h_t - mu_j)^2 / sigma2_j, up to
    // the constant that the normalisation removes anyway.
    arma::mat lq =
        arma::repmat(arma::log(mixture.weight) - 0.5 * arma::log(mixture.variance), tt, 1);
    lq -= 0.5 *
          arma::square(arma::repmat(y_star - h_col, 1, n) - arma::repmat(mixture.mean, tt, 1)) /
          arma::repmat(mixture.variance, tt, 1);

    arma::mat q = lq.each_col() - arma::max(lq, 1);
    q = arma::exp(q);
    q = q.each_col() / arma::sum(q, 1);

    // Inverse c.d.f. draw: the number of components whose cumulated probability
    // still exceeds the uniform draw counts down from the component that is hit.
    arma::uvec s =
        n - arma::sum(arma::repmat(arma::randu<arma::vec>(tt), 1, n) < arma::cumsum(q, 1), 1);

    // A uniform draw above the last cumulated probability -- reachable by
    // rounding alone, since the row sums to one only to within floating point --
    // leaves no component hit and would index one past the last. Round it down to
    // the last component rather than letting `.elem()` throw.
    s.elem(arma::find(s >= n)).fill(n - 1);

    return s;
}

/// Draws from \f$N(V^{-1} b, V^{-1})\f$ with V symmetric tridiagonal and
/// positive definite, given as its main diagonal and its constant off-diagonal.
///
/// V is a random walk precision plus a positive diagonal, so it carries one band
/// on either side of the main one and its Cholesky factor \f$V = R'R\f$ is upper
/// bidiagonal. Both come out of a single sweep over the periods, which is what
/// keeps the draw O(T). Assembling the T x T matrix and factorising it densely
/// -- the spelling this replaces -- costs O(T^2) memory and O(T^3) time for a
/// matrix that is zero everywhere but three diagonals, and at the sample lengths
/// a macroeconomic model is run on that is the bulk of the sampler.
///
/// One factorisation serves both the mean and the noise, and so does the back
/// substitution: \f$Rm = f\f$ and \f$Rn = z\f$ make \f$R(m + n) = f + z\f$, and
/// the sum is what the caller wants. With \f$Cov(R^{-1} z) = (R'R)^{-1}\f$ the
/// result has the posterior mean and the posterior covariance.
///
/// The factorisation is also the positive definiteness check. V is positive
/// definite by construction, being a precision plus a positive diagonal, so a
/// non-positive pivot means a conditioning value has already degenerated.
///
/// @param out T-vector, overwritten with the draw.
/// @param main_diag main diagonal of V.
/// @param off_diag the single off-diagonal value of V, the same on both sides.
/// @param rhs the vector b, so that the mean is \f$V^{-1} b\f$.
/// @param noise T-vector of independent standard normal draws.
/// @return `false` if V turns out not to be positive definite, in which case
///   `out` is meaningless.
inline bool draw_tridiagonal_normal(arma::vec &out, const arma::vec &main_diag, double off_diag,
                                    const arma::vec &rhs, const arma::vec &noise)
{
    const arma::uword tt = main_diag.n_elem;

    arma::vec r(tt);     // main diagonal of R
    arma::vec u(tt - 1); // the one super-diagonal of R
    arma::vec f(tt);     // R' f = b, by forward substitution

    double pivot = main_diag(0);
    if (!(pivot > 0.0))
    {
        return false;
    }
    r(0) = std::sqrt(pivot);
    f(0) = rhs(0) / r(0);

    for (arma::uword t = 1; t < tt; t++)
    {
        u(t - 1) = off_diag / r(t - 1);
        pivot = main_diag(t) - u(t - 1) * u(t - 1);
        if (!(pivot > 0.0))
        {
            return false;
        }
        r(t) = std::sqrt(pivot);
        f(t) = (rhs(t) - u(t - 1) * f(t - 1)) / r(t);
    }

    out.set_size(tt);
    out(tt - 1) = (f(tt - 1) + noise(tt - 1)) / r(tt - 1);
    for (arma::uword t = tt - 1; t-- > 0;)
    {
        out(t) = (f(t) + noise(t) - u(t) * out(t + 1)) / r(t);
    }

    return true;
}

} // namespace stochvol_detail

/**
 * @brief Draws the log-volatility of a stochastic volatility model.
 *
 * For a series of error terms \f$u_{it}\f$, \f$t = 1, \dots, T\f$, with
 * \f[
 *   u_{it} \sim N(0, \exp(h_{it})), \qquad
 *   h_{it} = h_{i,t-1} + v_{it}, \qquad v_{it} \sim N(0, \sigma_i),
 * \f]
 * the measurement equation is linearised by squaring and taking logarithms,
 * \f[
 *   \log(u_{it}^2 + c_i) = h_{it} + \varepsilon_{it},
 * \f]
 * where \f$\varepsilon_{it}\f$ is \f$\log \chi^2_1\f$ and \f$c_i\f$ is a small
 * offset that keeps the logarithm finite. That distribution is approximated by
 * `mixture`. Conditional on a mixture indicator per period the state space is
 * linear and Gaussian, so the whole path \f$h_{i1}, \dots, h_{iT}\f$ is drawn in
 * one block from its normal conditional posterior. Each column of `y` is handled
 * independently.
 *
 * The random walk is the specification both papers give, and it is deliberately
 * not parameterised. A unit transition and a fixed \f$h_{i0}\f$ are what make the
 * posterior precision tridiagonal with a *constant* off-diagonal, and what leave
 * the initial state in the first period alone -- the two facts the banded draw
 * below is built on. A state equation with a transition coefficient, a
 * time-varying innovation variance, or a random initial state is a different
 * model and belongs in a general state path sampler, not here: see
 * `kalman_durbin_koopman_2002`, which takes all three.
 *
 * @param algorithm the name the caller is known by, used to prefix a message.
 * @param mixture the normal mixture standing in for the log chi-squared error.
 * @param y T x K matrix of error terms, one period per row. Must be finite;
 *   zeros are permitted because of the offset in `constant`.
 * @param h T x K matrix of the current log-volatility, used as the conditioning
 *   value of the mixture indicators. Must have the dimensions of `y`.
 * @param sigma K-vector of the current variances of the log-volatility
 *   innovations. Must be finite and strictly positive.
 * @param h_init K-vector of the initial states \f$h_{i0}\f$. Must be finite.
 * @param constant K-vector of the offsets \f$c_i\f$ added before the logarithm.
 *   Must be finite and strictly positive.
 *
 * @return T x K matrix of the drawn log-volatility. The argument `h` is left
 *   untouched, so a caller may pass the same matrix it assigns the result to.
 *
 * @throws std::invalid_argument if the arguments do not describe a stochastic
 *   volatility model of at least two periods and one variable, if their
 *   dimensions disagree, or if any of them contains a non-finite or non-positive
 *   value where the algorithm needs a finite or positive one. The message names
 *   the argument.
 * @throws std::runtime_error if a posterior precision turns out not to be
 *   positive definite, or if a drawn path is not finite. Either means the
 *   conditioning values have degenerated, and the message says which variable it
 *   was.
 *
 * @warning The draw depends on the global Armadillo random number generator.
 *   Seed it with `arma::arma_rng::set_seed` for reproducible results.
 */
inline arma::mat stochvol_mixture_draw(const char *algorithm, const NormalMixture &mixture,
                                       const arma::mat &y, const arma::mat &h,
                                       const arma::vec &sigma, const arma::vec &h_init,
                                       const arma::vec &constant)
{
    using namespace stochvol_detail;

    // The tables cannot disagree if they were written as literals, but a mixture
    // assembled at runtime can, and the failure is a division by a zero variance
    // several steps later.
    require(algorithm, mixture.size() > 0, "the mixture has no components");
    require(algorithm, mixture.mean.n_elem == mixture.size() &&
                           mixture.variance.n_elem == mixture.size(),
            "the mixture's weights, means and variances differ in length");
    require(algorithm, mixture.variance.min() > 0.0,
            "the mixture has a component of non-positive variance");

    // Checks
    require(algorithm, y.n_rows == h.n_rows && y.n_cols == h.n_cols,
            "'h' must have the dimensions of 'y' (" + dims(y) + "), got " + dims(h));
    require(algorithm, y.n_cols > 0, "'y' must have at least one column");
    require(algorithm, y.n_rows > 1,
            "a stochastic volatility model needs at least two periods, got " +
                std::to_string(y.n_rows));
    require(algorithm, all_finite(y), "'y' contains NaN or infinite values");
    require(algorithm, all_finite(h), "'h' contains NaN or infinite values");

    const arma::uword k = y.n_cols;
    const arma::uword tt = y.n_rows;

    require_length(algorithm, sigma, k, "'sigma'");
    require_length(algorithm, h_init, k, "'h_init'");
    require_length(algorithm, constant, k, "'constant'");
    require_positive(algorithm, sigma, "'sigma'");
    require_positive(algorithm, constant, "'constant'");
    require(algorithm, all_finite(h_init), "'h_init' contains NaN or infinite values");

    // The random walk that the log-volatility follows enters the posterior as
    // the precision hh = D'D, with D the first difference operator: unit
    // diagonal, -1 below it. That product is tridiagonal -- 2 on the main
    // diagonal, 1 in the last period because no difference reaches past it, and
    // -1 on either side -- so only its two bands are formed and the T x T matrix
    // itself never is. Adding the mixture precisions leaves the posterior
    // precision tridiagonal as well, which is what the banded draw needs.
    arma::vec hh_diag(tt);
    hh_diag.fill(2.0);
    hh_diag(tt - 1) = 1.0;
    const double hh_off_diag = -1.0;

    arma::mat draw = h;

    for (arma::uword i = 0; i < k; i++)
    {
        const std::string variable = "variable " + std::to_string(i + 1);

        // Prepare series. Finite by construction: 'y' is finite and 'constant'
        // is strictly positive, so the argument of the logarithm is too.
        const arma::vec y_star = arma::log(arma::square(y.col(i)) + constant(i));

        // Sample the mixture indicators
        const arma::uvec s = draw_mixture_states(y_star, draw.col(i), mixture);

        // Sample the log-volatility in one block from N(V^-1 b, V^-1)
        const arma::vec mixture_precision = 1.0 / mixture.variance.elem(s);
        const arma::vec post_h_diag = hh_diag / sigma(i) + mixture_precision;
        const double post_h_off_diag = hh_off_diag / sigma(i);

        // hh 1 = e_1: every row of D'D sums to zero but the first, which sums to
        // one, so the initial state reaches the first period only.
        arma::vec rhs = mixture_precision % (y_star - mixture.mean.elem(s));
        rhs(0) += h_init(i) / sigma(i);

        arma::vec h_i;
        if (!draw_tridiagonal_normal(h_i, post_h_diag, post_h_off_diag, rhs,
                                     arma::randn<arma::vec>(tt)))
        {
            throw std::runtime_error(std::string(algorithm) +
                                     ": the posterior precision of the log-volatility of " +
                                     variable + " is not positive definite");
        }

        if (!all_finite(h_i))
        {
            throw std::runtime_error(std::string(algorithm) + ": the drawn log-volatility of " +
                                     variable + " is not finite");
        }

        draw.col(i) = h_i;
    }

    return draw;
}

} // namespace bayests::core

#endif // BAYESTS_CORE_ALGORITHMS_STOCHVOL_MIXTURE_H
