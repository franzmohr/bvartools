// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "core/algorithms/constrained_var_path.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * The data-completion step of a VAR whose panel is not observed whole.
 *
 * Given its parameters a VAR is a Gaussian over its own path, x = vec(Y), with
 * a banded precision: writing the model as H x = d + u, where H has identities
 * on its block diagonal and -A_j on the j-th block subdiagonal and d carries the
 * offsets and the presample, the precision is H' S^-1 H with S block diagonal,
 * and the canonical mean -- precision times mean -- is H' S^-1 d. No H is ever
 * formed: each period contributes one outer product over the p + 1 periods its
 * equation touches, which is the band.
 *
 * What was observed enters in three ways, and which way is decided by what the
 * row is, not by what frequency it describes:
 *
 *   - a hard row of one entry pins that entry. It is held rather than drawn:
 *     the precision is partitioned into the entries nothing pins and the rest,
 *     and the pinned part moves to the right hand side. When most of a panel is
 *     observed, most of the problem goes away here.
 *   - a soft row is an observation with an error, w'x + e = value, and adds
 *     w w' / variance to the precision like any Gaussian likelihood term. That
 *     can widen the band: a row over L periods couples entries L - 1 periods
 *     apart.
 *   - a hard row of several entries -- an aggregate observed exactly -- is not
 *     a likelihood term with a variance to add, so the draw is made without it
 *     and then conditioned on it (Rue and Held 2005, section 2.3.3): with V the
 *     precision's inverse applied to W', the draw z becomes
 *     z - V (W V)^-1 (W z - r). That holds W x = r to rounding in every draw
 *     and needs one extra pair of band solves per hard row.
 *
 * The log density is the other half and a different computation: a Kalman
 * filter over a state holding as many periods as the longer of p and the
 * widest row, which yields the density of each period's rows given everything
 * before them -- the pointwise form a model writes as its log likelihood. The
 * precision form would give the total at the same cost but not the terms.
 */

namespace bayests::core
{

namespace
{

void require(const bool ok, const std::string &what)
{
    if (!ok)
    {
        throw std::invalid_argument("constrained_var_path: " + what);
    }
}

std::string dims(const arma::mat &m)
{
    return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

/// True if every element is a finite number, read off the exponent bits, for
/// the reason chan_jeliazkov_2009.cpp gives: a host compiling with -ffast-math
/// may fold std::isfinite() to true.
bool all_finite(const arma::mat &m)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");
    for (arma::uword i = 0; i < m.n_elem; i++)
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

/// The sizes of the problem, read off the prior once.
struct Shape
{
    arma::uword k = 0;
    arma::uword tt = 0;
    arma::uword p = 0;
    bool a_stacked = false;
    bool s_stacked = false;
};

Shape check_prior(const VarPathPrior &prior)
{
    Shape s;
    s.k = prior.offset.n_rows;
    s.tt = prior.offset.n_cols;
    require(s.k > 0 && s.tt > 0,
            "'offset' must be k x T with both positive, got " + dims(prior.offset));
    require(all_finite(prior.offset), "'offset' must be finite");

    require(prior.coefficients.n_cols % s.k == 0,
            "'coefficients' must have a multiple of k = " + std::to_string(s.k) +
                " columns, one block per lag, got " + dims(prior.coefficients));
    s.p = prior.coefficients.n_cols / s.k;
    if (s.p > 0)
    {
        require(prior.coefficients.n_rows == s.k || prior.coefficients.n_rows == s.k * s.tt,
                "'coefficients' must have k = " + std::to_string(s.k) + " rows, or " +
                    std::to_string(s.k * s.tt) + " to give one block per period, got " +
                    dims(prior.coefficients));
        require(all_finite(prior.coefficients), "'coefficients' must be finite");
        s.a_stacked = prior.coefficients.n_rows != s.k;
        require(prior.presample.n_rows == s.k && prior.presample.n_cols == s.p,
                "'presample' must be k x p = " + std::to_string(s.k) + "x" +
                    std::to_string(s.p) + ", got " + dims(prior.presample));
        require(all_finite(prior.presample), "'presample' must be finite");
    }

    require(prior.covariance.n_cols == s.k &&
                (prior.covariance.n_rows == s.k || prior.covariance.n_rows == s.k * s.tt),
            "'covariance' must be k x k = " + std::to_string(s.k) + "x" + std::to_string(s.k) +
                ", or " + std::to_string(s.k * s.tt) + "x" + std::to_string(s.k) +
                " to give one block per period, got " + dims(prior.covariance));
    require(all_finite(prior.covariance), "'covariance' must be finite");
    s.s_stacked = prior.covariance.n_rows != s.k;
    return s;
}

/// Period t's block of an argument that is one block or a stack of them.
arma::mat period_block(const arma::mat &m, const arma::uword t, const arma::uword k,
                       const bool stacked)
{
    return stacked ? arma::mat(m.rows(t * k, t * k + k - 1)) : m;
}

arma::mat precision_of(const arma::mat &sigma)
{
    arma::mat out;
    if (!arma::inv_sympd(out, sigma))
    {
        throw std::invalid_argument(
            "constrained_var_path: 'covariance' is not symmetric positive definite");
    }
    return out;
}

/// The rows of a constraint set, regrouped: which entries each row has, and
/// the index each entry reaches in vec of the path.
struct Rows
{
    std::vector<std::vector<arma::uword>> entries;
    arma::uvec index;
};

Rows check_constraints(const Constraints &c, const Shape &s, const arma::vec &variances)
{
    const arma::uword n_rows = c.value.n_elem;
    const arma::uword n_entries = c.row.n_elem;

    require(c.group.n_elem == n_rows, "'group' must have one element per row");
    require(c.period.n_elem == n_entries && c.variable.n_elem == n_entries &&
                c.weight.n_elem == n_entries,
            "'row', 'period', 'variable' and 'weight' must have one element per entry");
    require(all_finite(c.value) && all_finite(c.weight), "values and weights must be finite");

    const arma::uword n_groups = n_rows > 0 ? c.group.max() : 0;
    require(variances.n_elem >= n_groups,
            "the rows name " + std::to_string(n_groups) + " soft groups, and " +
                std::to_string(variances.n_elem) + " variances were given");
    require(all_finite(variances) && (variances.n_elem == 0 || variances.min() > 0.0),
            "every soft variance must be positive and finite");

    Rows rows;
    rows.entries.assign(n_rows, {});
    rows.index.set_size(n_entries);
    for (arma::uword e = 0; e < n_entries; e++)
    {
        require(c.row(e) < n_rows, "an entry names a row past the last");
        require(c.period(e) < s.tt, "an entry names a period past the path");
        require(c.variable(e) < s.k, "an entry names a variable past the panel");
        rows.entries[c.row(e)].push_back(e);
        rows.index(e) = c.period(e) * s.k + c.variable(e);
    }
    for (arma::uword r = 0; r < n_rows; r++)
    {
        require(!rows.entries[r].empty(), "every row needs at least one entry");
    }
    return rows;
}

/// A symmetric matrix in lower band storage: data(i - j, j) is element (i, j)
/// for 0 <= i - j <= bw. After cholesky() the same storage holds the lower
/// factor.
struct Band
{
    arma::uword n = 0;
    arma::uword bw = 0;
    arma::mat data;

    Band(const arma::uword size, const arma::uword width)
        : n(size), bw(width), data(width + 1, size, arma::fill::zeros)
    {
    }

    /// Adds x to element (i, j) and, the matrix being symmetric, to (j, i).
    void add(arma::uword i, arma::uword j, const double x)
    {
        if (i < j)
        {
            std::swap(i, j);
        }
        data(i - j, j) += x;
    }
};

/// Overwrites a band with its lower Cholesky factor. False if a pivot is not
/// positive, which is the matrix not being positive definite.
bool cholesky(Band &b)
{
    for (arma::uword j = 0; j < b.n; j++)
    {
        const arma::uword k0 = j > b.bw ? j - b.bw : 0;
        double d = b.data(0, j);
        for (arma::uword q = k0; q < j; q++)
        {
            const double l = b.data(j - q, q);
            d -= l * l;
        }
        if (!(d > 0.0))
        {
            return false;
        }
        d = std::sqrt(d);
        b.data(0, j) = d;

        const arma::uword last = std::min(b.n - 1, j + b.bw);
        for (arma::uword i = j + 1; i <= last; i++)
        {
            // Unsigned: i - bw is only taken where it cannot wrap round.
            const arma::uword k1 = std::max(k0, i > b.bw ? i - b.bw : arma::uword(0));
            double sum = b.data(i - j, j);
            for (arma::uword q = k1; q < j; q++)
            {
                sum -= b.data(i - q, q) * b.data(j - q, q);
            }
            b.data(i - j, j) = sum / d;
        }
    }
    return true;
}

/// Solves L y = x in place.
void solve_lower(const Band &l, arma::vec &x)
{
    for (arma::uword i = 0; i < l.n; i++)
    {
        const arma::uword k0 = i > l.bw ? i - l.bw : 0;
        double sum = x(i);
        for (arma::uword q = k0; q < i; q++)
        {
            sum -= l.data(i - q, q) * x(q);
        }
        x(i) = sum / l.data(0, i);
    }
}

/// Solves L' y = x in place.
void solve_upper(const Band &l, arma::vec &x)
{
    for (arma::uword i = l.n; i-- > 0;)
    {
        const arma::uword last = std::min(l.n - 1, i + l.bw);
        double sum = x(i);
        for (arma::uword q = i + 1; q <= last; q++)
        {
            sum -= l.data(q - i, i) * x(q);
        }
        x(i) = sum / l.data(0, i);
    }
}

/// Everything the draw and the moments share: the factor of the precision of
/// the entries nothing pins, the mean before the hard rows, and what
/// conditioning on the hard rows needs.
struct Conditioned
{
    arma::uword n = 0;
    std::vector<bool> pinned;
    arma::vec held;        // the pinned values, by index in vec of the path
    arma::uvec free_index; // index in vec of the path of each free entry
    Band factor{0, 0};
    arma::vec mean;        // of the free entries, before the hard rows
    arma::mat w;           // hard rows of several entries, over the free entries
    arma::vec r;           // what they have to equal, the pinned part moved over
    arma::mat v;           // precision^-1 w'
    arma::mat c_root;      // upper Cholesky factor of w v
};

/// The width of the band: the VAR's own, p + 1 periods of k entries, or the
/// widest soft row if that is wider. Hard rows of several entries never enter
/// the precision and so do not count.
arma::uword band_width(const Constraints &c, const Shape &s, const Rows &rows)
{
    arma::uword bw = (s.p + 1) * s.k - 1;
    for (arma::uword r = 0; r < c.value.n_elem; r++)
    {
        if (c.group(r) == 0)
        {
            continue;
        }
        arma::uword lo = s.k * s.tt, hi = 0;
        for (const arma::uword e : rows.entries[r])
        {
            lo = std::min(lo, rows.index(e));
            hi = std::max(hi, rows.index(e));
        }
        bw = std::max(bw, hi - lo);
    }
    return bw;
}

/// The VAR's precision over the whole path and its canonical mean, H' S^-1 H
/// and H' S^-1 d, added into `precision` and `canonical`.
void add_var_prior(const VarPathPrior &prior, const Shape &s, Band &precision,
                   arma::vec &canonical)
{
    const arma::uword k = s.k;
    const arma::uword p = s.p;
    const arma::mat constant_precision =
        s.s_stacked ? arma::mat() : precision_of(prior.covariance);
    const arma::mat identity = arma::eye<arma::mat>(k, k);

    for (arma::uword t = 0; t < s.tt; t++)
    {
        const arma::mat a = p > 0 ? period_block(prior.coefficients, t, k, s.a_stacked) : arma::mat();
        const arma::mat s_inv =
            s.s_stacked ? precision_of(period_block(prior.covariance, t, k, true))
                        : constant_precision;

        // The part of period t's mean that is not a lag inside the path.
        arma::vec d = prior.offset.col(t);
        for (arma::uword j = 1; j <= p; j++)
        {
            if (t < j)
            {
                d += a.cols((j - 1) * k, j * k - 1) * prior.presample.col(p + t - j);
            }
        }

        // g[j] is the coefficient of y_{t-j} in u_t: the identity for j = 0,
        // -A_j after it, and only as far back as the path goes.
        const arma::uword reach = std::min<arma::uword>(p, t);
        std::vector<arma::mat> g(reach + 1);
        g[0] = identity;
        for (arma::uword j = 1; j <= reach; j++)
        {
            g[j] = -a.cols((j - 1) * k, j * k - 1);
        }

        const arma::vec s_inv_d = s_inv * d;
        for (arma::uword j = 0; j <= reach; j++)
        {
            canonical.subvec((t - j) * k, (t - j) * k + k - 1) += g[j].t() * s_inv_d;
        }
        for (arma::uword j = 0; j <= reach; j++)
        {
            const arma::mat s_inv_g = s_inv * g[j];
            for (arma::uword i = 0; i <= j; i++)
            {
                // Rows in period t - i, columns in period t - j: on or below the
                // diagonal, and for i == j the one triangle of a symmetric block.
                const arma::mat block = g[i].t() * s_inv_g;
                for (arma::uword col = 0; col < k; col++)
                {
                    for (arma::uword row = (i == j ? col : 0); row < k; row++)
                    {
                        precision.add((t - i) * k + row, (t - j) * k + col, block(row, col));
                    }
                }
            }
        }
    }
}

Conditioned condition(const VarPathPrior &prior, const Constraints &c,
                      const arma::vec &variances)
{
    const Shape s = check_prior(prior);
    const Rows rows = check_constraints(c, s, variances);
    const arma::uword n = s.k * s.tt;
    const arma::uword n_rows = c.value.n_elem;

    Band precision(n, band_width(c, s, rows));
    arma::vec canonical(n, arma::fill::zeros);
    add_var_prior(prior, s, precision, canonical);

    // Soft rows, as likelihood terms on the whole path.
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (c.group(r) == 0)
        {
            continue;
        }
        const double variance = variances(c.group(r) - 1);
        const std::vector<arma::uword> &entries = rows.entries[r];
        for (arma::uword a = 0; a < entries.size(); a++)
        {
            const arma::uword ea = entries[a];
            canonical(rows.index(ea)) += c.weight(ea) * c.value(r) / variance;
            for (arma::uword b = 0; b <= a; b++)
            {
                const arma::uword eb = entries[b];
                // Each unordered pair once into the one stored triangle; a pair
                // of distinct entries reaching the same element counts twice.
                const double scale = (a != b && rows.index(ea) == rows.index(eb)) ? 2.0 : 1.0;
                precision.add(rows.index(ea), rows.index(eb),
                              scale * c.weight(ea) * c.weight(eb) / variance);
            }
        }
    }

    // What the hard rows of one entry pin.
    Conditioned out;
    out.n = n;
    out.pinned.assign(n, false);
    out.held.zeros(n);
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (c.group(r) == 0 && rows.entries[r].size() == 1)
        {
            const arma::uword e = rows.entries[r][0];
            require(c.weight(e) != 0.0, "a hard row of one entry needs a weight other than zero");
            out.pinned[rows.index(e)] = true;
            out.held(rows.index(e)) = c.value(r) / c.weight(e);
        }
    }

    std::vector<arma::uword> position(n, 0);
    std::vector<arma::uword> free_index;
    for (arma::uword i = 0; i < n; i++)
    {
        if (!out.pinned[i])
        {
            position[i] = free_index.size();
            free_index.push_back(i);
        }
    }
    out.free_index = arma::uvec(free_index);
    const arma::uword n_free = free_index.size();

    // The precision of the free entries, and their canonical mean with the
    // pinned ones moved to the right hand side.
    Band reduced(n_free, precision.bw);
    arma::vec rhs(n_free);
    for (arma::uword f = 0; f < n_free; f++)
    {
        rhs(f) = canonical(out.free_index(f));
    }
    for (arma::uword j = 0; j < n; j++)
    {
        const arma::uword last = std::min(n - 1, j + precision.bw);
        for (arma::uword i = j; i <= last; i++)
        {
            const double value = precision.data(i - j, j);
            if (value == 0.0)
            {
                continue;
            }
            const bool free_i = !out.pinned[i];
            const bool free_j = !out.pinned[j];
            if (free_i && free_j)
            {
                reduced.data(position[i] - position[j], position[j]) += value;
            }
            else if (free_i)
            {
                rhs(position[i]) -= value * out.held(j);
            }
            else if (free_j)
            {
                rhs(position[j]) -= value * out.held(i);
            }
        }
    }

    // Hard rows of several entries, over the free entries.
    std::vector<arma::uword> hard;
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (c.group(r) == 0 && rows.entries[r].size() > 1)
        {
            hard.push_back(r);
        }
    }
    out.w.zeros(hard.size(), n_free);
    out.r.zeros(hard.size());
    for (arma::uword h = 0; h < hard.size(); h++)
    {
        out.r(h) = c.value(hard[h]);
        for (const arma::uword e : rows.entries[hard[h]])
        {
            const arma::uword i = rows.index(e);
            if (out.pinned[i])
            {
                out.r(h) -= c.weight(e) * out.held(i);
            }
            else
            {
                out.w(h, position[i]) += c.weight(e);
            }
        }
    }

    if (n_free == 0)
    {
        return out;
    }

    if (!cholesky(reduced))
    {
        throw std::runtime_error("constrained_var_path: the precision of the entries nothing "
                                 "pins is not positive definite");
    }
    out.factor = reduced;
    out.mean = rhs;
    solve_lower(out.factor, out.mean);
    solve_upper(out.factor, out.mean);

    if (!hard.empty())
    {
        out.v = out.w.t();
        for (arma::uword h = 0; h < hard.size(); h++)
        {
            arma::vec column = out.v.col(h);
            solve_lower(out.factor, column);
            solve_upper(out.factor, column);
            out.v.col(h) = column;
        }
        arma::mat wv = out.w * out.v;
        wv = 0.5 * (wv + wv.t());
        if (!arma::chol(out.c_root, wv))
        {
            throw std::runtime_error(
                "constrained_var_path: the hard rows of several entries are not independent "
                "given the rest -- two of them determine the same combination");
        }
    }
    return out;
}

/// (W V)^-1 x through the Cholesky factor of W V.
arma::vec solve_hard(const Conditioned &c, const arma::vec &x)
{
    const arma::vec half = arma::solve(arma::trimatl(c.c_root.t()), x);
    return arma::solve(arma::trimatu(c.c_root), half);
}

/// The free entries conditioned on the hard rows: z - V (W V)^-1 (W z - r).
arma::vec impose_hard(const Conditioned &c, const arma::vec &z)
{
    if (c.w.n_rows == 0)
    {
        return z;
    }
    return z - c.v * solve_hard(c, c.w * z - c.r);
}

arma::mat assemble(const Conditioned &c, const arma::vec &free, const arma::uword k)
{
    arma::vec x = c.held;
    for (arma::uword f = 0; f < c.free_index.n_elem; f++)
    {
        x(c.free_index(f)) = free(f);
    }
    return arma::reshape(x, k, c.n / k);
}

} // namespace

arma::mat draw_constrained_var_path(const VarPathPrior &prior, const Constraints &constraints,
                                    const arma::vec &soft_variances)
{
    const Conditioned c = condition(prior, constraints, soft_variances);
    const arma::uword k = prior.offset.n_rows;
    if (c.free_index.n_elem == 0)
    {
        return assemble(c, arma::vec(), k);
    }

    arma::vec noise = arma::randn<arma::vec>(c.free_index.n_elem);
    solve_upper(c.factor, noise);
    return assemble(c, impose_hard(c, c.mean + noise), k);
}

ConstrainedPathMoments constrained_var_path_moments(const VarPathPrior &prior,
                                                    const Constraints &constraints,
                                                    const arma::vec &soft_variances)
{
    const Conditioned c = condition(prior, constraints, soft_variances);
    const arma::uword k = prior.offset.n_rows;
    const arma::uword n_free = c.free_index.n_elem;

    ConstrainedPathMoments out;
    out.covariance.zeros(c.n, c.n);
    if (n_free == 0)
    {
        out.mean = assemble(c, arma::vec(), k);
        return out;
    }
    out.mean = assemble(c, impose_hard(c, c.mean), k);

    arma::mat inverse = arma::eye<arma::mat>(n_free, n_free);
    for (arma::uword f = 0; f < n_free; f++)
    {
        arma::vec column = inverse.col(f);
        solve_lower(c.factor, column);
        solve_upper(c.factor, column);
        inverse.col(f) = column;
    }
    if (c.w.n_rows > 0)
    {
        const arma::mat root_inv = arma::solve(arma::trimatl(c.c_root.t()), c.v.t());
        inverse -= root_inv.t() * root_inv;
    }
    out.covariance.submat(c.free_index, c.free_index) = inverse;
    return out;
}

arma::vec constrained_var_path_log_density(const VarPathPrior &prior,
                                           const Constraints &constraints,
                                           const arma::vec &soft_variances)
{
    const Shape s = check_prior(prior);
    const Rows rows = check_constraints(constraints, s, soft_variances);
    const arma::uword k = s.k;
    const arma::uword p = s.p;
    const arma::uword n_rows = constraints.value.n_elem;

    // Where each row ends, and how many periods the state has to carry for the
    // widest of them to be a function of it.
    std::vector<std::vector<arma::uword>> ending(s.tt);
    arma::uword span = std::max<arma::uword>(p, 1);
    for (arma::uword r = 0; r < n_rows; r++)
    {
        arma::uword lo = s.tt, hi = 0;
        for (const arma::uword e : rows.entries[r])
        {
            lo = std::min(lo, constraints.period(e));
            hi = std::max(hi, constraints.period(e));
        }
        ending[hi].push_back(r);
        span = std::max(span, hi - lo + 1);
    }
    const arma::uword m = k * span;

    // The state before the path: its p known periods, and zeros beyond them,
    // which nothing reads -- no row reaches before the path and the transition
    // only p periods back.
    arma::vec state(m, arma::fill::zeros);
    for (arma::uword j = 1; j <= std::min(p, span); j++)
    {
        state.subvec((j - 1) * k, j * k - 1) = prior.presample.col(p - j);
    }
    arma::mat state_cov(m, m, arma::fill::zeros);

    arma::mat transition(m, m, arma::fill::zeros);
    if (span > 1)
    {
        transition.submat(k, 0, m - 1, m - k - 1) = arma::eye<arma::mat>(m - k, m - k);
    }

    arma::vec out(s.tt, arma::fill::zeros);
    const double log_two_pi = std::log(2.0 * arma::datum::pi);

    for (arma::uword t = 0; t < s.tt; t++)
    {
        if (p > 0)
        {
            transition.submat(0, 0, k - 1, k * p - 1) =
                period_block(prior.coefficients, t, k, s.a_stacked);
        }
        const arma::mat sigma = period_block(prior.covariance, t, k, s.s_stacked);

        state = transition * state;
        state.head(k) += prior.offset.col(t);
        state_cov = transition * state_cov * transition.t();
        state_cov.submat(0, 0, k - 1, k - 1) += sigma;

        const std::vector<arma::uword> &now = ending[t];
        if (now.empty())
        {
            continue;
        }

        const arma::uword count = now.size();
        arma::mat z(count, m, arma::fill::zeros);
        arma::vec observed(count), noise(count, arma::fill::zeros);
        for (arma::uword i = 0; i < count; i++)
        {
            const arma::uword r = now[i];
            observed(i) = constraints.value(r);
            if (constraints.group(r) > 0)
            {
                noise(i) = soft_variances(constraints.group(r) - 1);
            }
            for (const arma::uword e : rows.entries[r])
            {
                const arma::uword lag = t - constraints.period(e);
                z(i, lag * k + constraints.variable(e)) += constraints.weight(e);
            }
        }

        const arma::vec innovation = observed - z * state;
        const arma::mat pz = state_cov * z.t();
        arma::mat f = z * pz;
        f.diag() += noise;
        f = 0.5 * (f + f.t());

        arma::mat root;
        if (!arma::chol(root, f))
        {
            throw std::runtime_error(
                "constrained_var_path: what period " + std::to_string(t + 1) +
                " observes is determined by what was observed before it, or twice over");
        }
        const arma::vec scaled = arma::solve(arma::trimatl(root.t()), innovation);
        out(t) = -0.5 * (static_cast<double>(count) * log_two_pi +
                         2.0 * arma::accu(arma::log(root.diag())) + arma::dot(scaled, scaled));

        // gain' = F^-1 Z P, through the factor
        const arma::mat gain_t =
            arma::solve(arma::trimatu(root), arma::solve(arma::trimatl(root.t()), pz.t()));
        state += gain_t.t() * innovation;
        state_cov -= pz * gain_t;
        state_cov = 0.5 * (state_cov + state_cov.t());
    }
    return out;
}

} // namespace bayests::core
