# Vendored BayesTS core

`src/core/` and `inst/include/bayests/` are a copy of the core layer of
[BayesTS](https://github.com/franzmohr/BayesTS): the sampler declarations in
`include/bayests/` and the numerics in `src/core/`. Nothing else of that
project is here -- the core deliberately links neither HDF5 nor HighFive,
prints nothing and reads no files, which is what makes it embeddable in an R
package at all.

The copy is **BayesTS `4a64082`**, upstream `main` after the `v0.3.0` release
(tagged on `18a86c2`). Past the release it carries the fix to
`core/models/var_tvp_discount.cpp` (`2dc9250`), which is what makes
`add_predictive_loglik()` score more than the first horizon of a discounted VAR;
the flat-prior warning for `bvs` (`fa6dd89`, `f2abd4e`), which adds
`flat_selection_prior()` and `flat_selection_message()` to `priors.h` and reaches
R through `RcppReporter::warnings()` and `.raise_core_warnings()`; the three
SSVS refusals in `validate()` (`a85abe2`); and the non-centred random walks of
`VarTvpStochvol` (`59c495f`), which bring `core/models/noncentred_support.h`
and reach R through `omega_v` in `add_priors()` and the `omega*` draws
`src/VarTvpStochvol.cpp` returns. The previous refresh sat at `6fe91d2`, 0.3.0
plus the discount fix. 0.3.0 is not archived yet, so
there is no version DOI to name; the concept DOI
<https://doi.org/10.5281/zenodo.22722531> resolves to the newest release
whenever one is deposited, and the last archived release is 0.2.0,
<https://doi.org/10.5281/zenodo.22765348>.

Upstream commits that change nothing under `include/` or `src/core/` do not move
the copy off that commit; a refresh that copies anything newer has to update this
paragraph. Nothing enforces that -- the script compares the vendored *files* and
`inst/COPYRIGHTS`, never this prose -- and the same paragraph in dfmtools, which
vendors part of the same core, spent two refreshes naming a version it was no
longer at. So check this one against the upstream `git log` on the way out, as
part of the refresh rather than after it.

Upstream layout is preserved, so a refresh is a copy of two directories plus
the patch below -- which upstream has since made unnecessary, though the script
still applies it. `tools/update-bayests-core.R` does both:

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

Keep this list short; every entry is something a refresh has to reapply. It is
down to one, and that one no longer needs reapplying.

1. **Armadillo comes from RcppArmadillo** -- which is now upstream's own
   arrangement rather than something reapplied here. Every copied header and
   source reaches Armadillo through `bayests/arma.h`, and that header includes
   whatever `BAYESTS_ARMA_HEADER` names, defaulting to `<armadillo>` for a
   standalone build. `src/Makevars` and `src/Makevars.win` define it as
   `<RcppArmadillo.h>`, which is what points Armadillo's RNG at R's, so
   `set.seed()` reaches the samplers.

   It began as a patch. `arma.h` existed only here, and the refresh script
   rewrote every

   ```cpp
   #define ARMA_DONT_PRINT_FAST_MATH_WARNING
   #include <armadillo>
   ```

   into an include of it. Upstream has since adopted both the header -- byte for
   byte the same file -- and the rule, so nothing under `include/bayests/` or
   `src/core/` includes `<armadillo>` any more and the rewrite has nothing left
   to do. It is kept all the same, as is the script's refusal to finish while a
   copied file still includes `<armadillo>`: neither costs anything, and either
   would catch a file upstream forgot. The script's `keep` list still names
   `arma.h`, which is now merely redundant -- upstream's copy is identical, so
   copying it over ours is a no-op.

   What is load-bearing is the define. Lose it and the vendored core breaks
   silently: it compiles, links and runs, and merely stops honouring
   `set.seed()`. `src/bayests_rng_guard.cpp` fails the build with a
   `static_assert` if Armadillo ends up configured for its own RNG, which is the
   one thing standing between a dropped `-D` and a package whose seeds do
   nothing.

That is the whole list. `core/models/var_tvp_wishart.cpp` used to be a second
entry -- it was not copied, because upstream neither listed it in
`src/core/CMakeLists.txt` nor defined the `VarTvpWishartInput::validate()` its
sampler calls. Both are fixed upstream and the file is vendored like any other.

## Not copied

`skip` in the refresh script holds the upstream sources this package does not
take. Each one needs the reason written here.

**The dynamic factor models.** Upstream has four of them --
`DfmNormalGamma`, `DfmNormalStochvol`, `DfmTvpGamma` and `DfmTvpStochvol` -- and
all four are dfmtools' models rather than this package's. Their headers and
sources are skipped, along with `core/models/dfm_support.h`, which is the shared
half none of them compiles without. Nothing here includes any of them: the R side
of the DFM is gone and there was never a `src/Dfm*.cpp` binding, so they compiled
into the shared object with nothing able to reach them.

Dropping them is safe to link precisely because nothing includes them: the only
`Dfm` symbols any other translation unit defines are the four
`Dfm*Input::validate()` methods in `core/inputs.cpp`, and each of those throws on
bad input and calls nothing from the skipped files.

**A new DFM upstream has to be added to `skip` by hand, and forgetting does more
than bloat the object.** `dfm_support.h` is skipped, so a sampler copied without
it does not compile at all -- which is how the omission surfaces, at the next
build rather than at the refresh. The refresh script says how many DFMs upstream
has beside the list, so the two can be compared.

**The factor augmented VAR.** `FavarNormalWishart` is a factor model, so it is
dfmtools' rather than this package's, and it is skipped along with
`core/models/favar_support.h`. The mechanism is not merely the same as the DFMs',
it is the DFMs': `favar_support.h` includes the `dfm_support.h` skipped just
above, so a copy of the sampler without them would fail to compile rather than
merely sit unreachable in the shared object. Upstream has one of these.

Everything the paragraph below says about the DFMs' type surface holds for it
too. `FavarNormalWishartInitial` and `FavarNormalWishartInput` are in
`bayests/inputs.h`, `FavarNormalWishartDraws` in `bayests/results.h`,
`n_obs_factors`, `n_state()`, `n_favar_lambda()` and `n_favar_a()` in
`bayests/spec.h`, `TrainData::f_obs` in `bayests/data.h`, and
`FavarNormalWishartInput::validate()` in `core/inputs.cpp` -- all files every
model shares and which are copied whole. That validator calls only the local
helpers in `inputs.cpp`, so it links with the sampler absent, exactly as the four
DFM validators do.

**The factor model's Kalman filter.** `core/models/factor_score.h` is what
scores a forecast of a factor model: its history reaches the density through the
latent factors, so the realised observation of each scored period updates their
distribution before the next is predicted, and the score is a prediction error
decomposition rather than a log likelihood over another sample. Only the four
DFM sources include it, so here it would be unreachable.

It is the one skipped file that would compile if it were copied: unlike the
samplers above it reaches for nothing but `bayests/spec.h` and
`core/models/predictive_score.h`, both of which are here. It is skipped all the
same, because `inst/COPYRIGHTS` has to name every vendored file and a reviewer
reading that list should not find one that nothing includes. The VAR and VEC
models are scored by `predictive_score.h` alone -- with realised history their
regressors do not depend on the draw, so the score is the model's own pointwise
log likelihood over the scored periods, and no filter is needed.

`core/algorithms/chan_jeliazkov_2009.cpp` is copied and now carries a second
entry point, `chan_jeliazkov_2009_conditional`, which holds the trailing elements
of every state column at observed values instead of drawing them. Nothing here
calls it -- it exists for the FAVAR -- and it costs one function in a file this
package already vendors for the time-varying models. Upstream verified the split
that introduced it left every fingerprint unmoved.

This does **not** take the DFMs out of the vendored core, and cannot. Every
`Dfm*Initial` and `Dfm*Input` is in `bayests/inputs.h`, every `Dfm*Draws` in
`bayests/results.h`, `n_factors` and `n_lambda()` in `bayests/spec.h`, and the
four validators in `core/inputs.cpp`. Those four files are shared by every model
and are copied whole, so the DFMs' type surface and their input validation stay
compiled in. What goes is the samplers that would act on them.
`src/bayests_r_io.h` reads `n_factors` for the same reason: the field is there
whether or not anything sets it.

## The discounted models

`core/models/var_tvp_discount.cpp`, `core/models/vec_tvp_discount.cpp` and the
`core/models/discount_support.h` both of them share arrived with BayesTS 0.3.0
and were left out of the previous refresh, under a *Not copied* entry that said
to delete itself once something here could reach them. `src/VarTvpDiscount.cpp`
and `src/VecTvpDiscount.cpp` are that something, so they are vendored like any
other model and the entry is gone.

**They are the one pair in the core that is not a sampler.** The posterior is
the matrix normal dynamic linear model of West & Harrison (1997, ch. 16) with
the discounted Wishart of Uhlig (1997), closed form in one pass over the sample,
which is why the bindings look different from the eight beside them:

* The entry point is `estimate()` returning a `Var/VecTvpDiscountPosterior`,
  not `draw_coefficients()` returning draws. What crosses into R is one column
  per period -- `a$mean`, `a$scale`, `a$cov`, `u_sigma$scale` and `df` -- and it
  carries no `mcpar`, because there is no chain to have a start, an end or a
  thinning interval. `draws_to_r()` still does the transposing, the orientation
  being the same one a chain crosses in.
* Estimation consumes no random numbers, so a model estimated here and the same
  model estimated by the `bayests` command line agree to the bit rather than
  merely in distribution. Forecasting does consume them, being i.i.d. draws from
  the closed form, and runs under the model's seed like every sampler.
* The log-likelihood entry point runs the filter again rather than reading the
  stored posterior back. `log_likelihood()` returns `posterior.loglik`, which
  `estimate()` fills on its way through and which is not among the blocks the
  posterior is written as. Re-running costs one deterministic pass over the
  sample and is the same arithmetic on the same input, not a second estimate of
  it -- the same reasoning upstream's own file front-end gives.

`VarSpec::delta_beta` and `VarSpec::delta_sigma` had to be added to
`read_spec()` in `src/bayests_r_io.h`. Both default to one in the core, and one
is a model rather than a neutral value -- the quantity the discount governs does
not move -- so a specification whose discounts never reached the filter would
have estimated a constant coefficient model and said nothing about it.

The matrix normal prior is read by `read_matrix_normal_prior()` from `mean` and
`cov`, deliberately not from the `mu` and `v_inv` a sampler's normal prior over
the vectorised coefficients uses. The two are different objects: `cov` is the
regressor side of a covariance in units of the error covariance, the full prior
covariance being `Sigma` kronecker `cov`, which is what makes the posterior
closed form.

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

All eight VAR algorithms are converted: `VarNormalWishart`, `VarNormalGamma`,
`VarNormalStochvol`, `VarNormalAld`, `VarTvpGamma`, `VarTvpWishart`,
`VarTvpStochvol` and `VarTvpAld`. Nothing in `src/` samples any more; the
numerics all live in `src/core/`.

The seven VEC algorithms are converted too -- `VecNormalWishart`,
`VecNormalGamma`, `VecNormalStochvol`, `VecKlgs2010`, `VecTvpGamma`,
`VecTvpStochvol` and `VecTvpWishart` -- so `src/bvecalg.cpp`, this package's own
implementation of the VEC, is gone. `VecNormalWishart` was the last of them to
wait on a binding, and while it did the two were separate implementations of the
same model with only the vendored one carrying the fixes recorded in BayesTS.
One implementation each now.

Nothing here is now vendored without being reachable. `VecNormalWishart` held
that position until it got a binding, `VarNormalAld` and `VarTvpAld` until
theirs, and the DFM samplers until they moved to `skip`; see *Not copied*.

The two quantile VARs are the exception to the file-per-model rule above in one
respect: each has `.Var*AldCoefficients` and `.Var*AldLogLik` and no
`.Var*AldForecasts`, because the sampler's `forecast()` always throws. Three
things about them a reader of those two files should know, all of them silent
when wrong:

1. **The quantile has to cross.** `read_spec()` reads `model$quantile` into
   `VarSpec::quantile`, which defaults to 0.5. Lose that line and a model asked
   for the 0.8 quantile estimates the median, runs to completion and says
   nothing -- which is the failure `test-quantile_var.R` checks the estimand
   against rather than the shape of the output.
2. **The refusals are the estimand's, not gaps.** No covariance block, because
   rotating the errors makes the estimand a combination of quantiles rather than
   the quantile of a combination; no forecast, because the `h` step quantile is
   not the quantile of the iterated one step quantiles. `validate()` rejects
   both, the bindings let it, and `add_posterior_forecasts()` refuses before
   reaching C++ so the message names the model rather than the horizon.
3. **`u_sigma_inv` is read for its column count alone.** The log likelihood of
   these models is the asymmetric Laplace density, which is closed form and
   marginal of the latent scales, so nothing reads the precision path for its
   values -- but `iterations()` counts its columns. The two loglik readers take
   the last period only, which is `k * k` numbers per draw instead of
   `k * k * tt`.

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
