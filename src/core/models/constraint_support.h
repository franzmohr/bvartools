// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_CORE_MODELS_CONSTRAINT_SUPPORT_H
#define BAYESTS_CORE_MODELS_CONSTRAINT_SUPPORT_H

#include "bayests/data.h"
#include "bayests/spec.h"
#include "core/models/model_support.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace bayests::core
{

/// How far an observed value may be from the entry of `y` it pins and still
/// count as the same number, relative to the larger of one and that entry.
///
/// A host writes both from one source, so the two agree to the bit when the
/// weight is one and to a rounding error when it is not. What this has to stop
/// is a row written against the wrong period or the wrong variable, and those
/// are off by the size of the data.
constexpr double kPinnedValueTolerance = 1e-10;

/// Refuses a constraint set that does not describe a panel of k variables over
/// tt periods, or that could not be conditioned on.
///
/// `where` names the group in the file, for the messages. Positions in them are
/// one-based, as the file counts them.
///
/// Beyond the shapes, four rules, each of which a sampler would otherwise
/// meet as something other than a refusal:
///
///   - every row carries at least one entry, and one weight that is not zero.
///     A row of zeros observes nothing and, if hard, states 0 = value;
///   - no (row, period, variable) appears twice, since whether the two weights
///     were meant to be added is the file's to say, not the reader's;
///   - soft groups are numbered 1..G with none left out, so that group g is
///     the g-th variance a sampler draws and a gap is not a variance with no
///     rows to inform it;
///   - a hard row of one entry pins that entry, and it has to agree with what
///     `y` holds there -- a row written against the wrong period is off by the
///     size of the data. Two hard rows may not pin one entry, and a hard row of
///     several entries has to reach at least one entry nothing pins: otherwise
///     it is either implied or contradicted by the others, and conditioning on
///     it exactly is a singular system.
///
/// `held` is the panel the pins have to agree with, k x tt, or null where there
/// is none -- a scenario for the horizon pins values nothing has observed -- and
/// then only the shape k x tt is checked against.
inline void validate_constraints(const Constraints &c, const arma::uword k, const arma::uword tt,
                                 const arma::mat *held, const std::string &where)
{
    const arma::uword n_rows = c.value.n_elem;
    const arma::uword n_entries = c.row.n_elem;

    if (c.group.n_elem != n_rows)
    {
        throw std::invalid_argument(where + "/group must have one element per row of " + where +
                                    "/value, " + std::to_string(n_rows) + ", got " +
                                    std::to_string(c.group.n_elem));
    }
    if (c.period.n_elem != n_entries || c.variable.n_elem != n_entries ||
        c.weight.n_elem != n_entries)
    {
        throw std::invalid_argument(
            where + "/row, period, variable and weight must have one element per entry each, got " +
            std::to_string(n_entries) + ", " + std::to_string(c.period.n_elem) + ", " +
            std::to_string(c.variable.n_elem) + " and " + std::to_string(c.weight.n_elem));
    }
    require_finite(c.value, where + "/value");
    require_finite(c.weight, where + "/weight");

    if (k == 0 || tt == 0)
    {
        throw std::invalid_argument(where + " constrains a panel, and there is none: the "
                                            "observations it belongs to are empty");
    }

    std::vector<arma::uword> entries_per_row(n_rows, 0);
    std::vector<bool> nonzero_in_row(n_rows, false);
    for (arma::uword e = 0; e < n_entries; e++)
    {
        if (c.row(e) >= n_rows)
        {
            throw std::invalid_argument(where + "/row holds " + std::to_string(c.row(e) + 1) +
                                        ", but there are " + std::to_string(n_rows) + " rows");
        }
        if (c.period(e) >= tt)
        {
            throw std::invalid_argument(where + "/period holds " + std::to_string(c.period(e) + 1) +
                                        ", but the panel has " + std::to_string(tt) + " periods");
        }
        if (c.variable(e) >= k)
        {
            throw std::invalid_argument(where + "/variable holds " +
                                        std::to_string(c.variable(e) + 1) + ", but the panel has " +
                                        std::to_string(k) + " variables");
        }
        entries_per_row[c.row(e)]++;
        if (c.weight(e) != 0.0)
        {
            nonzero_in_row[c.row(e)] = true;
        }
    }
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (entries_per_row[r] == 0 || !nonzero_in_row[r])
        {
            throw std::invalid_argument(where + " row " + std::to_string(r + 1) +
                                        " has no entry with a weight other than zero, so it "
                                        "observes nothing");
        }
    }

    // Duplicates, found by sorting the triplets rather than by a set: the
    // sizes are those of a panel, and a sort needs nothing beside the index.
    const arma::uword panel = k * tt;
    std::vector<arma::uword> order(n_entries);
    for (arma::uword e = 0; e < n_entries; e++)
    {
        order[e] = e;
    }
    const auto key = [&](const arma::uword e) {
        return c.row(e) * panel + c.period(e) * k + c.variable(e);
    };
    std::sort(order.begin(), order.end(),
              [&](const arma::uword a, const arma::uword b) { return key(a) < key(b); });
    for (arma::uword i = 1; i < n_entries; i++)
    {
        if (key(order[i]) == key(order[i - 1]))
        {
            const arma::uword e = order[i];
            throw std::invalid_argument(where + " row " + std::to_string(c.row(e) + 1) +
                                        " names period " + std::to_string(c.period(e) + 1) +
                                        ", variable " + std::to_string(c.variable(e) + 1) +
                                        " twice");
        }
    }

    // Soft groups, 1..G with none left out.
    if (n_rows > 0)
    {
        const arma::uword n_groups = c.group.max();
        std::vector<bool> seen(n_groups + 1, false);
        for (arma::uword r = 0; r < n_rows; r++)
        {
            seen[c.group(r)] = true;
        }
        for (arma::uword g = 1; g <= n_groups; g++)
        {
            if (!seen[g])
            {
                throw std::invalid_argument(
                    where + "/group runs to " + std::to_string(n_groups) + " but holds no row of group " +
                    std::to_string(g) + "; soft groups are numbered from one with none left out");
            }
        }
    }

    // What the hard single-entry rows pin, and that it is what y holds there.
    std::vector<bool> pinned(panel, false);
    std::vector<arma::uword> first_entry(n_rows, n_entries);
    for (arma::uword e = 0; e < n_entries; e++)
    {
        first_entry[c.row(e)] = std::min(first_entry[c.row(e)], e);
    }
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (c.group(r) != 0 || entries_per_row[r] != 1)
        {
            continue;
        }
        const arma::uword e = first_entry[r];
        const arma::uword t = c.period(e);
        const arma::uword i = c.variable(e);
        if (pinned[t * k + i])
        {
            throw std::invalid_argument(where + " pins period " + std::to_string(t + 1) +
                                        ", variable " + std::to_string(i + 1) +
                                        " in two hard rows; one observation is one row");
        }
        pinned[t * k + i] = true;

        if (held == nullptr)
        {
            continue;
        }
        const double observed = c.value(r) / c.weight(e);
        const double value = (*held)(i, t);
        if (!(std::abs(observed - value) <= kPinnedValueTolerance * std::max(1.0, std::abs(value))))
        {
            throw std::invalid_argument(
                where + " row " + std::to_string(r + 1) + " observes " + std::to_string(observed) +
                " at period " + std::to_string(t + 1) + ", variable " + std::to_string(i + 1) +
                ", where the observations hold " + std::to_string(value) +
                "; a row written against the wrong period or variable reads as this");
        }
    }
    std::vector<bool> reaches_free(n_rows, false);
    for (arma::uword e = 0; e < n_entries; e++)
    {
        if (!pinned[c.period(e) * k + c.variable(e)])
        {
            reaches_free[c.row(e)] = true;
        }
    }
    for (arma::uword r = 0; r < n_rows; r++)
    {
        if (c.group(r) != 0 || entries_per_row[r] < 2)
        {
            continue;
        }
        if (!reaches_free[r])
        {
            throw std::invalid_argument(
                where + " row " + std::to_string(r + 1) +
                " is hard, and every entry it reaches is pinned by another hard row, so it is "
                "either implied by them or contradicts them");
        }
    }
}

/// validate_constraints() against a panel that is there, `by_period` being
/// k x tt, one column per period: its shape, and what the pins have to agree
/// with.
inline void validate_constraints(const Constraints &c, const arma::mat &by_period,
                                 const std::string &where)
{
    validate_constraints(c, by_period.n_rows, by_period.n_cols, &by_period, where);
}

/// Refuses constraints on an algorithm that does not read them, and checks
/// them on one that does.
///
/// Refused rather than ignored: the entries a panel does not observe hold
/// placeholders, and a sampler that reads `y` whole would estimate a model from
/// them and report it as a fit to data. `supported` is the algorithm's own
/// answer and `algorithm` names it, as in require_supported_iid_block().
///
/// The training constraints are checked against `y` stacked by period, which
/// is the one layout its three spellings share (stacked_response()); the test
/// constraints against `/data/test/y`, one row per period.
inline void require_supported_constraints(const VarSpec &spec, const TrainData &train,
                                          const TestData &test, const bool supported,
                                          const char *algorithm)
{
    if (train.constraints.empty() && test.constraints.empty())
    {
        return;
    }

    if (!supported)
    {
        throw std::invalid_argument(
            std::string(algorithm) +
            " does not read /data/train/constraints or /data/test/constraints: it estimates from a "
            "panel observed whole, and would read the values standing in for what was not "
            "observed as data. Drop the constraints, or fill the panel");
    }

    const arma::uword k = spec.k > 0 ? static_cast<arma::uword>(spec.k) : 0;
    if (!train.constraints.empty())
    {
        const arma::uword tt = k > 0 ? train.y.n_elem / k : 0;
        const arma::mat by_period =
            tt > 0 ? arma::mat(arma::reshape(stacked_response(train), k, tt)) : arma::mat();
        validate_constraints(train.constraints, by_period, "/data/train/constraints");
    }
    if (!test.constraints.empty())
    {
        validate_constraints(test.constraints, arma::trans(test.y), "/data/test/constraints");
    }
}

/// The same for a scenario, /data/forecast/constraints: refused where the
/// algorithm does not condition its forecast, and otherwise checked against a
/// horizon of `h` periods. Nothing has been observed there, so the pins are
/// checked for their shape and not against any value.
inline void require_supported_forecast_constraints(const VarSpec &spec,
                                                   const ForecastData &forecast,
                                                   const bool supported, const char *algorithm)
{
    if (forecast.constraints.empty())
    {
        return;
    }

    if (!supported)
    {
        throw std::invalid_argument(
            std::string(algorithm) +
            " does not read /data/forecast/constraints: it does not condition a forecast on a "
            "scenario. Drop the constraints for an unconditional forecast");
    }
    if (spec.h <= 0)
    {
        throw std::invalid_argument(
            "/data/forecast/constraints conditions a forecast, and /model asks for none: set h");
    }

    validate_constraints(forecast.constraints, static_cast<arma::uword>(spec.k),
                         static_cast<arma::uword>(spec.h), nullptr, "/data/forecast/constraints");
}

} // namespace bayests::core

#endif // BAYESTS_CORE_MODELS_CONSTRAINT_SUPPORT_H
