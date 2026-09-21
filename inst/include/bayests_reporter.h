// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef BVARTOOLS_BAYESTS_REPORTER_H
#define BVARTOOLS_BAYESTS_REPORTER_H

#include "bayests/arma.h"
#include "bayests/reporter.h"

#include <string>
#include <vector>

namespace bvartools
{

/// The R half of the BayesTS `Reporter` contract.
///
/// The core prints nothing and cancels nothing on its own -- it cannot, because
/// a package may neither write to stdout nor call R's API from wherever it
/// pleases. Everything it wants to say arrives here instead, and this is the
/// only place in the R build that talks to the console.
///
/// The samplers call `check_interrupt()` and `progress()` once per draw and do
/// no throttling of their own, so both are throttled here: `Rcpp::checkUserInterrupt()`
/// is an R evaluation and a chain of 20000 draws should not pay for 20000 of
/// them. `interrupt_every` matches the every-hundredth-draw check the
/// hand-written samplers in this package already use.
class RcppReporter final : public bayests::Reporter
{
public:
    explicit RcppReporter(bool verbose = false, long long interrupt_every = 100)
        : verbose_(verbose), interrupt_every_(interrupt_every > 0 ? interrupt_every : 1)
    {
    }

    /// A line the core marks `"warning: "` is kept, verbose or not, for the
    /// binding to hand back to R once the sampler has returned -- see
    /// `warnings()`. Everything else is status and printed only when verbose.
    ///
    /// Not raised here with `Rcpp::warning()`: that is `Rf_warning()`, which
    /// longjmps out under `options(warn = 2)` or a `tryCatch(warning = )`, and a
    /// longjmp out of the middle of a sampler skips the destructors of every
    /// Armadillo matrix it holds.
    void message(const std::string &text) override
    {
        static const std::string prefix = "warning: ";
        if (text.compare(0, prefix.size(), prefix) == 0)
        {
            warnings_.push_back(text.substr(prefix.size()));
            return;
        }
        if (verbose_)
        {
            Rcpp::Rcout << text << std::endl;
        }
    }

    /// The warnings the run raised, without their prefix, in the order it
    /// raised them. A binding returns these as the `warnings` element of its
    /// result and `.raise_core_warnings()` turns each into an R warning.
    const std::vector<std::string> &warnings() const { return warnings_; }

    /// Reports at most once per percent, and only when asked to be verbose.
    void progress(long long done, long long total) override
    {
        if (!verbose_ || total <= 0)
        {
            return;
        }

        const int percent = static_cast<int>((done * 100) / total);
        if (percent == last_percent_ && done != total)
        {
            return;
        }
        last_percent_ = percent;

        Rcpp::Rcout << "\rProgress: " << percent << "% (" << done << "/" << total << ")";
        if (done == total)
        {
            Rcpp::Rcout << std::endl;
        }
    }

    void finish() override
    {
        if (verbose_ && last_percent_ >= 0)
        {
            Rcpp::Rcout << std::endl;
            last_percent_ = -1;
        }
    }

    /// Throws `Rcpp::internal::InterruptedException` when the user has pressed
    /// Ctrl-C, which unwinds the sampler and is turned back into an R condition
    /// by the `BEGIN_RCPP` / `END_RCPP` block around the binding.
    void check_interrupt() override
    {
        if (++calls_ % interrupt_every_ == 0)
        {
            Rcpp::checkUserInterrupt();
        }
    }

private:
    bool verbose_;
    long long interrupt_every_;
    long long calls_ = 0;
    int last_percent_ = -1;
    std::vector<std::string> warnings_;
};

} // namespace bvartools

#endif // BVARTOOLS_BAYESTS_REPORTER_H
