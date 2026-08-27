# Vendored BayesTS core

`src/core/` and `inst/include/bayests/` are a copy of the core layer of
[BayesTS](https://github.com/franzmohr/BayesTS): the sampler declarations in
`include/bayests/` and the numerics in `src/core/`. Nothing else of that
project is here -- the core deliberately links neither HDF5 nor HighFive,
prints nothing and reads no files, which is what makes it embeddable in an R
package at all.

Upstream layout is preserved, so a refresh is a copy of two directories plus
the patch below. `tools/update-bayests-core.R` does both:

```bash
Rscript tools/update-bayests-core.R path/to/BayesTS/source_tree
```

It reports what it changed, what it skipped and anything here that upstream no
longer has, and it is idempotent -- running it twice changes nothing the second
time. Everything below is what that script encodes; keep the two in step.

| Upstream | Here |
| --- | --- |
| `include/bayests/*.h` | `inst/include/bayests/*.h` |
| `src/core/{spec,inputs}.cpp` | `src/core/` |
| `src/core/models/*` | `src/core/models/` |
| `src/core/algorithms/*` | `src/core/algorithms/` |

`core/algorithms/kalman_durbin_koopman_2002.cpp` and
`core/algorithms/stochvol_ocsn_2007.cpp` replaced the package's own copies of
those functions; see *The simulation smoother* and *The stochastic volatility
draw* below.

The copied files keep their `SPDX-License-Identifier: BSD-3-Clause` headers.
That is compatible with this package's GPL (>= 2) and the copyright holder is
the same person, so there is nothing to reconcile. `inst/COPYRIGHTS` spells the
difference out for a CRAN reviewer and lists the files by name. That list has to
match what is actually vendored, and the refresh script compares the two and
refuses to finish when they disagree -- it does not edit the list, so a refresh
that added or removed a file stops until the list is brought into line by hand.

## Local modifications

Keep this list short; every entry is something a refresh has to reapply.

1. **Armadillo comes from RcppArmadillo.** Upstream opens with

   ```cpp
   #define ARMA_DONT_PRINT_FAST_MATH_WARNING
   #include <armadillo>
   ```

   and every such `#include <armadillo>` is replaced by
   `#include "bayests/arma.h"`, a header that exists only here and includes
   `<RcppArmadillo.h>`. That is what points Armadillo's RNG at R's, so
   `set.seed()` reaches the samplers. The rule is uniform -- every copied file
   that includes Armadillo, header or source -- which is why it can be applied
   by a script rather than kept as a list of file names that goes stale.

   Losing this patch is the only way to break the vendored core silently: it
   compiles, links and runs, and merely stops honouring `set.seed()`. Two things
   watch for it. The script refuses to finish while any copied file still
   includes `<armadillo>`, and `src/bayests_rng_guard.cpp` fails the build with
   a `static_assert` if Armadillo ends up configured for its own RNG.

That is the whole list. `core/models/var_tvp_wishart.cpp` used to be a second
entry -- it was not copied, because upstream neither listed it in
`src/core/CMakeLists.txt` nor defined the `VarTvpWishartInput::validate()` its
sampler calls. Both are fixed upstream and the file is vendored like any other.

## Not copied

`skip` in the refresh script holds the upstream sources this package does not
take. Each one needs the reason written here.

**The dynamic factor model.** `bayests/dfm_normal_gamma.h`,
`core/models/dfm_normal_gamma.cpp` and `core/models/dfm_support.h` are dfmtools'
model, not this package's. Nothing here includes any of the three -- the R side
of the DFM is gone and there was never a `src/DfmNormalGamma.cpp` binding -- so
they compiled into the shared object with nothing able to reach them. 636 lines
of the 12,592 vendored.

Dropping them is safe to link precisely because nothing includes them: the only
`Dfm` symbol any other translation unit defines is
`DfmNormalGammaInput::validate()`, in `core/inputs.cpp`, and that function
throws on bad input and calls nothing from the three files.

This does **not** take the DFM out of the vendored core, and cannot.
`DfmNormalGammaInitial` and `DfmNormalGammaInput` are in `bayests/inputs.h`,
`DfmNormalGammaDraws` in `bayests/results.h`, `n_factors` and `n_lambda()` in
`bayests/spec.h`, and the validator above in `core/inputs.cpp`. Those four are
shared by every model and are copied whole, so the DFM's type surface and its
input validation stay compiled in. What goes is the sampler that would act on
them. `src/bayests_r_io.h` reads `n_factors` for the same reason: the field is
there whether or not anything sets it.

## The simulation smoother

`core/algorithms/kalman_durbin_koopman_2002.cpp` replaces the package's own
copy of the same function, which used to live in
`src/kalman_durbin_koopman_2002.cpp` and was character-for-character identical.
That file is gone, and with it the `bvartools::kalman_durbin_koopman_2002`
entry point that `// [[Rcpp::interfaces(cpp)]]` used to generate in
`inst/include/bvartools_RcppExports.h` — a package that linked against it has
to call the core function instead. It was never exposed to R, so nothing in
`NAMESPACE` or `man/` changed at the time.

The R side has followed since. `src/kalman_dk.cpp` held a *third* copy of the
smoother — the same numerics again, line for line, differing only in taking its
arguments by value where the core takes four by reference — and was the exported
R function `kalman_dk()`. It is gone. `kalman_durbin_koopman_2002()` replaces it,
a thin `Rcpp::export` wrapper over the core function in
`src/kalman_durbin_koopman_2002.cpp`, so R and the samplers run the same code.
Because the two bodies were identical, a draw from a given seed is bit-identical
to what `kalman_dk()` returned; nothing about any posterior changes.

The wrapper carries no `// [[Rcpp::interfaces(r, cpp)]]`, so `bvartools::kalman_dk`
is gone from `inst/include/bvartools_RcppExports.h` too. That entry point was
added deliberately — *"Made `kalman_dk` callable from C++"*, NEWS 0.2.2 — so this
is a removal, not a tidy-up, and an external package that used it has to link the
core source instead. The reason is the one this file spends a section on: a
`bvartools::` entry point wraps every call in an `Rcpp::RNGScope`, and a caller
that holds the RNG across calls gets the rewind demonstrated below. Regenerating
it under the name `kalman_durbin_koopman_2002` would have been worse still, since
that is a name this package has already retired for exactly that reason.

Its documented example was rewritten on the way. It called `gen_var()`, which no
longer exists, and reached for `temp$data$Y` and `temp$data$SUR`, which the model
object stopped carrying under those names some time before that.

`VarTvpGammaCoefficients.cpp`, `VarTvpStochvolCoefficients.cpp` and
`VarTvpWishartCoefficients.cpp` now call it as a plain C++ function. They used
to call `bvartools::kalman_durbin_koopman_2002`, i.e. this package's own C++
interface, which reaches the same code through `R_GetCCallable` after wrapping
every argument in a `SEXP` — once per draw.

**This changes the draws, and the new ones are the correct ones.** The
generated interface wraps each call in an `Rcpp::RNGScope`, which calls
`GetRNGstate()` on entry and `PutRNGstate()` on exit. Inside a sampler that
already holds the RNG, that entry re-reads `.Random.seed` and *discards
everything drawn since the previous interface call*:

```
        three uniforms drawn in a row, with one interface call in between
        no call:    0.2655087  0.3721239  0.5728534
        with call:  0.2655087  0.2655087  0.3721239   <- the stream rewound
```

So in every iteration the smoother rewound the stream to where the previous
smoother call had left it, and the numbers the intervening `a0` draw had
consumed were handed out a second time. The first iteration of a chain is
unaffected — it is bit-identical before and after this change — and from the
second one the two disagree. The smoother itself is not in question: called on
the same inputs from the same seed, the vendored source returns bit-identical
values and leaves the RNG in the same state as the old route.

The package had two more call sites through its own interface, and they carried
the same rewind: `bvartools::sur_const_to_tvp` in `post_bvs.cpp` and
`bvartools::stochvol_ocsn2007_internal` in `bvecalg.cpp`. Both now call the
function directly, the first declared in `src/sur_const_to_tvp.h` and the second
in `core/algorithms/stochvol_ocsn_2007.h`. A deterministic callee was no
protection — the rewind demonstrated above used one — so `post_bvs` drew
differently too.

**Nothing inside this package may call a `bvartools::` function.** Those
wrappers exist for other packages; from in here they are the same code with a
rewound RNG. `grep -rn 'bvartools::' src/*.cpp` should stay empty.

## The stochastic volatility draw

`core/algorithms/stochvol_ocsn_2007.cpp` replaces `src/stochvol_ocsn2007_internal.cpp`,
which held the same algorithm. That file is gone, and with it the
`bvartools::stochvol_ocsn2007_internal` entry point that `// [[Rcpp::interfaces(cpp)]]`
used to generate in `inst/include/bvartools_RcppExports.h` — a package that
linked against it has to call the core function instead. It was never exposed to
R, so nothing in `NAMESPACE` or `man/` changed at the time.

It is exposed to R now. `stochvol_ocsn_2007()` and `stochvol_ksc_1998()` are
thin `Rcpp::export` wrappers over the two core functions, in
`src/stochvol_ocsn_2007.cpp` and `src/stochvol_ksc_1998.cpp`, and they replace
the exported R implementations `stochvol_ocsn2007()` and `stochvol_ksc1998()`,
which are gone. Those were separate implementations of the same two algorithms
and had to be kept in step by hand — which is how `stochvol_ksc1998()` came to
be missing the two fixes its sibling received. One implementation each now, in
the core, and R reaches it rather than reimplementing it.

Two things about those wrappers. Their C++ functions are named
`stochvol_*_export` and exported under the core function's name via
`// [[Rcpp::export(stochvol_ocsn_2007)]]`, because the wrapper and the function
it calls cannot both be `stochvol_ocsn_2007` in C++. And they carry no
`// [[Rcpp::interfaces(r, cpp)]]`: a `bvartools::` entry point is what rewound
the RNG in the three cases above, and nothing needs one here — a package that
links against this one can call the core function.

`R/algo_bvectvp.R` called `stochvol_ksc1998()` and now calls
`stochvol_ksc_1998()`. Its draws from a given seed change; that sampler is not
reachable through `add_priors()` yet, which refuses a VEC TVP prior, so the
change is not observable from R at present.

The vendored routine differs from the copy it replaces in two ways beyond the
name.

It reports bad input by throwing, which `END_RCPP` turns into an R error, and
the checks cover what the old copy indexed on trust: the lengths of `sigma`,
`h_init` and `constant` against the columns of `y`, and each of them being
finite, with `sigma` and `constant` strictly positive. `sigma(i) == 0` used to
divide by zero and `-DARMA_NO_DEBUG` in `src/Makevars` meant a short `sigma`
read past the end of its vector rather than raising anything.

It draws the mixture indicator in logs, shifted by the row maximum. Weighting
`arma::normpdf` and dividing by the row sum — what the old copy did — underflows
to ten zeros once an observation sits far enough out in the tails of all ten
components; the row then normalises to NaN and the indicator comes out as one
past the last component. With `-DARMA_NO_DEBUG`, `.elem()` used to read past the
end of the mixture instead of throwing. The R function had the same defect and
returned a matrix of `NA` for it.

Both spellings pick the same component whenever the densities do not underflow:
over 200,000 periods driven by the same uniforms they agreed exactly. Draws of
the VEC sampler with stochastic volatility still move by a rounding error,
because the posterior mean now comes from the same Cholesky factor as the draw
instead of from a separate LU solve — the same substitution
`core/models/model_support.h` documents upstream.

Since then the draw of the path has been rewritten upstream to use a banded
Cholesky. The posterior precision is tridiagonal — the random walk contributes
the two bands of `D'D`, the mixture only a diagonal — and the earlier spelling
materialised it as a dense `T x T` matrix and factorised that, at `O(T^3)` time
and `O(T^2)` memory. The factor of a symmetric tridiagonal is bidiagonal, so it
now comes out of one sweep over the periods at `O(T)`. Upstream measured 4.6x at
`T = 200` and 49x at `T = 1000` against the dense path, and verified the two
agree to a relative 1.8e-15 from the same seed, so draws are unchanged bar
rounding. Both algorithms share the routine, in
`core/algorithms/stochvol_mixture.h`, and each of the two `.cpp` files is its
mixture table plus a call; upstream's `test/unit_stochvol.cpp` covers them.

`core/algorithms/stochvol_ksc_1998.cpp` arrived with that refresh and is the
seven-component mixture of Kim, Shephard and Chib (1998). Nothing in this
package's samplers uses it — `bvecalg.cpp` and `core/models/var_tvp_stochvol.cpp`
both take the ten-component one — but `stochvol_ksc_1998()` exposes it to R.

## Building

`src/Makevars` adds `-I. -I../inst/include`, sets `CXX_STD = CXX17` and extends
`OBJECTS` with `core/*.cpp core/*/*.cpp`, since R only picks up sources at the
top level of `src/`. That last part uses `$(wildcard)`, hence
`SystemRequirements: GNU make` in `DESCRIPTION`. Nothing has to be listed by
name: a source added to a refresh is picked up by the next build.

## Using it from R

The core is values in, values out, and knows nothing about `SEXP`. Translation
is a layer of its own, the R counterpart of BayesTS's `src/io/hdf5/`, which is
not vendored because it translates from files rather than from R:

| | |
| --- | --- |
| `src/bayests_r_io.h` | the parts every model shares: spec, priors, positions, the transposes |
| `src/<Model>.cpp` | one file per model: its readers and writer in an anonymous namespace, then the three `[[Rcpp::export]]` entry points |

All six VAR algorithms are converted: `VarNormalWishart`, `VarNormalGamma`,
`VarNormalStochvol`, `VarTvpGamma`, `VarTvpWishart` and `VarTvpStochvol`.
Nothing in `src/` samples any more; the numerics all live in `src/core/`.

`VecNormalWishart` is the exception, and now the only thing here that is
vendored without being reachable -- the DFM sampler was the other one until it
moved to `skip`; see *Not copied*. The difference is that this one is waiting on
a binding rather than belonging to another package. Its sampler,
`core/algorithms/vec_to_var.*` and
`bayests/vec_normal_wishart.h` are all copied and compiled into the shared
object, but there is no `src/VecNormalWishart.cpp` binding and no R entry point,
so the VEC still runs through this package's own `src/bvecalg.cpp`. Until a
binding exists the two are separate implementations of the same model, and only
the vendored one has the fixes recorded in BayesTS -- among them the
cointegration space prior's scalar shrinkage, the loadings being written back
after reparameterisation, and the lag ordering of a simulated forecast path.

The R pipeline above all of this did not change: `bvarpost()`,
`add_posterior_forecasts()` and `add_posterior_loglik()` still dispatch on
`model$algorithm` to the same `.<Model>{Coefficients,Forecasts,LogLik}` names,
which still take and return the same lists.

Three things a binding has to get right, all of them silent when wrong:

1. **Draws are transposed at the boundary.** The samplers accumulate draws down
   the columns; R keeps them in rows. Everything crossing goes through
   `draws_to_r()` or `read_draws_if_present()`.
2. **Positions are one-based in R and zero-based in the core.**
   `read_positions()` is the only place that converts, and it rejects a zero
   rather than letting it wrap.
3. **`train$y` is already stacked** as `vec(y')`. The core reads a single
   column the same way, so it goes across untouched.

Pass a `bvartools::RcppReporter` (`inst/include/bayests_reporter.h`) for
progress and a throttled `Rcpp::checkUserInterrupt()`. Let the core's
exceptions out: `validate()` names the first inconsistency it finds and Rcpp
turns that into an R error, which beats the dimension mismatch a bad input used
to surface as.
