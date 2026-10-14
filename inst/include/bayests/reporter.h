// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#ifndef BAYESTS_REPORTER_H
#define BAYESTS_REPORTER_H

#include <string>

#include "bayests/arma.h"

namespace bayests
{

/// Everything a sampler wants to tell the outside world, and the one thing it
/// needs to ask.
///
/// The samplers run for minutes and have to stay usable from a command line, a
/// notebook and an R session, none of which agree on where text goes. Writing
/// to stdout directly settles that question in the numeric code, which is both
/// wrong for an embedded caller and, for an R package, not permitted -- R owns
/// its own console and a package that bypasses it is rejected.
///
/// check_interrupt() is the part that cannot be emulated afterwards. A host
/// that wants a long chain to be cancellable needs a call inside the draw loop
/// to cancel from, and it signals cancellation the only way it can from a
/// callback: by throwing.
class Reporter
{
public:
    virtual ~Reporter() = default;

    /// A line of human-readable status. Not a diagnostic: failures are thrown.
    virtual void message(const std::string &text) { (void)text; }

    /// Called once per draw. Implementations are expected to throttle -- the
    /// samplers do not, because how often a host can afford to redraw is the
    /// host's business.
    virtual void progress(long long done, long long total) { (void)done; (void)total; }

    /// The run has ended, successfully or not. Somewhere to close a progress
    /// bar or emit the trailing newline.
    virtual void finish() {}

    /// Throws to abort the run. Called on the same schedule as progress().
    virtual void check_interrupt() {}

    /// A kept draw of a block, as soon as it is made.
    ///
    /// For a host that shows the posterior forming while the chain runs -- the
    /// bands of a time-varying coefficient narrowing as draws accumulate --
    /// rather than once it has been written. `block` names the output group
    /// (`"a"` for /posterior/a), `index` is the draw's zero-based column in it,
    /// and `value` is the block in the sampler's own layout: for a time-varying
    /// block one row per coefficient and one column per period. The reference
    /// is valid for the duration of the call only; copy what is kept.
    ///
    /// Called only for kept draws, after burn-in and thinning, and so far only
    /// by the samplers with a time-varying coefficient path. Consumes no random
    /// numbers, and the default does nothing.
    virtual void kept_draw(const std::string &block, long long index, const arma::mat &value)
    {
        (void)block;
        (void)index;
        (void)value;
    }
};

/// Discards everything and never interrupts. The default for callers that just
/// want the draws, and for tests.
class NullReporter final : public Reporter
{
};

} // namespace bayests

#endif // BAYESTS_REPORTER_H
