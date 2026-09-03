# bvartools (development version)

* **`spillover()`**, the connectedness measures of Diebold and Yilmaz (2012):
  the total spillover index, the directional spillovers to and from each
  variable, their net difference and the net pairwise table. Methods for
  `bvarmodel`, `modellist` and `expandingwindow`, the last of which gives the
  index over a growing sample -- the chart that literature is built around --
  from the windows `use_expanding_window()` already produces. VEC models go
  through `vec_to_var()` first, as they do for `fevd()`.

    Every measure is computed once per posterior draw and only then summarised,
    which is not the same as computing it from the mean decomposition table: the
    measures are ratios, so the index of the mean is not the mean of the index.
    It also gives the index a credible interval, which the point estimate based
    original does not have. `keep_draws = TRUE` returns the draws themselves.

    The generalised decomposition it uses is the one of Pesaran and Shin (1998),
    which divides by the variance of the *shock* variable. This is **not** what
    `fevd(type = "gir")` computes -- see the next entry -- so `spillover()`
    carries its own worker, `.spillover_table()`, rather than calling
    `.vardecomp()`. That worker also fills the whole k x k table in one pass,
    where looping `.vardecomp()` over responses would rebuild the impulse
    responses once per row.

* **Fixed: `fevd()` read the wrong period of a time-varying or stochastic
  volatility model.** **Draws are unchanged; the variance decomposition of those
  models changes, and the old numbers were wrong.** The sample length was taken
  as `nrow(train$y) / k`, which was right when `gen_var()` stacked the series
  into one column and is wrong for the one-row-per-period layout
  `create_bvarmodel()` produces. With `k` variables it understated the sample by
  a factor of `k`, so the default `period` was fractional -- indexing the
  posterior at a truncated row -- and any `period` past `tt / k` was refused as
  "implausible" although it was in the sample. Constant-coefficient models are
  unaffected: they never read `tt`.

    The count now comes from the regressor matrix, whose row count is
    unambiguously `k` times the sample length, and falls back on `y` only where
    a model has no regressors, handling both layouts. Shared by `fevd()` and
    `spillover()` through a new internal `.collect_draws()`, lifted out of
    `fevd.bvarmodel()` so the two cannot disagree about which slice of a row a
    period is.

* **Known inconsistency, not changed.** `fevd(type = "gir")` divides by
  `sqrt(sigma_jj)` of the *response* variable (`src/vardecomp.cpp`), where
  Pesaran and Shin (1998) divide by `sigma_kk` of the *shock*. The package's own
  documentation states a third thing, `1 / sigma_jj` of the response without the
  square root. The factor cancels when a row is normalised, so
  `normalise_gir = TRUE` output is unaffected by the discrepancy between the
  documented and the coded version -- but neither equals Pesaran and Shin, and
  the two coincide only when every variable has unit residual variance.
  `bgvars::gfevd()` carries the same line. Correcting it would move existing
  `fevd()` and `gfevd()` numbers, so it is left for a deliberate decision;
  `test-spillover.R` pins the current behaviour in both directions so that a
  change to it fails loudly.

* **Vendored BayesTS core refreshed.** **Draws are unchanged**, for every VAR and
  VEC model this package samples. Upstream's own fingerprint comparison was run
  over the change and reports 74 fixtures unchanged and none moved, twice: once
  over the one piece of shared code it touched and once over the sampler added on
  top of it.

    Nothing that arrived is reachable from here. Upstream added two dynamic
    factor models, `DfmTvpGamma` and `DfmTvpStochvol`, which belong to
    `dfmtools`; both are in the refresh script's `skip` list beside the two that
    were already there, so neither is copied. What does arrive is what the two
    added to the files every model shares -- their `Input`, `Initial` and `Draws`
    structs in `bayests/inputs.h` and `bayests/results.h`, and their
    `validate()` methods in `core/inputs.cpp` -- which is the same arrangement
    the first two DFMs have had since they moved to `skip`, and which compiles
    their type surface in without their samplers.

    Beyond that: `core/models/model_support.h` gains `draw_random_walk_state()`,
    a helper the dynamic factor models share and nothing here calls, and three
    files in `core/algorithms/` gain a corrected citation in a comment -- Kim,
    Shephard and Chib rather than Kohn. No line of numerics changed.

    One thing worth knowing for the next refresh: the skip list now has to be
    extended by hand whenever upstream adds a dynamic factor model, and
    forgetting is not merely wasteful. `core/models/dfm_support.h` is skipped, so
    a sampler copied without it fails to compile, which is how this refresh found
    the two new ones. `src/core/VENDORED.md` says so beside the list.

# bvartools 1.0.0

* Added a `testthat` test suite, run by `R CMD check` and so by the existing
  R-CMD-check workflow. 591 assertions over 25 files, about ten seconds end to
  end, covering the workflow from `create_bvarmodel` and `create_bvecmodel`
  through priors, initial values, posterior simulation, model selection,
  impulse responses, forecast error variance decomposition, forecasts and the
  HDF5 round trip, plus the exported numerical helpers on their own. Fitted
  models are cached across files in `tests/testthat/helper-fixtures.R`, so the
  samplers run once per session rather than once per file.

  The tests check arithmetic wherever there is arithmetic to check, rather than
  only that a call returns without error: `loglik_normal` against the
  multivariate normal density written out, `minnesota_prior` against the
  Minnesota variance formula term by term, the information criteria against
  their definitions, the impulse response against the `Phi_h = A_1 Phi_{h-1}`
  recursion run by hand over the draws, `vec_to_var` against
  `A_1 = alpha beta' + I`, the forecast error variance decomposition against the
  shares summing to one, `post_normal` and the two variance samplers against the
  analytic posterior they draw from, and `covar_vector_to_matrix`,
  `generate_lower_block_diagonal`, `sur_const_to_tvp` and
  `coint_prepare_sur_data` against their block structures element by element.

  Four tests are written but skipped, each naming the defect that keeps it
  skipped, so that fixing the defect is a matter of removing one line.
  `add_initial_values(method = "prior")` reads the error variance prior from
  `priors$sigma`, which `add_priors` has never written -- it writes
  `priors$u_sigma` -- and `k` and `tt` are bound only in the `"ols"`/`"maxlik"`
  branch of `.add_initial_values_measurement_errors`, so the branch cannot run
  for any error specification or either model class. The same helper derives
  `tt` from the stacked `KT` vector of endogenous variables, so a stochastic
  volatility model gets an initial log volatility of `KT x K` and the sampler
  rejects it, which `add_posterior_coefficients` reports by returning the model
  with `error = TRUE`; `error = "sv"` therefore cannot currently be estimated
  through this workflow. `predict.bvarmodel` clamps `n_ahead` against
  `model$h` but then builds its result from `model$h` throughout, so a shorter
  horizon than was simulated is silently ignored. And
  `read_expanding_window_model_from_folder` calls `import_model_from_hdf5`,
  which no longer exists, and labels its result `expandwindmodellist` rather
  than the `expandingwindow` class `write_to_hdf5` took, so no
  `expandingwindow` method would dispatch on it even once the call is fixed.

  Two smaller things the tests record as they stand rather than skip. The
  guard in `combine_models` cannot fire: `classes[-which(classes == "list")]`
  drops every element when no input carries a `"list"` class, so an
  unsupported object passes the check and fails further down on
  `class(result) <- ...` of `NULL`. And `Matrix` is not declared in
  `DESCRIPTION` although `covar_vector_to_matrix`,
  `generate_lower_block_diagonal`, `sur_const_to_tvp`, `covar_prepare_data`,
  `coint_prepare_sur_data`, `post_bvs` and the two `post_gamma_*` functions all
  return or require its sparse classes; it is added to `Suggests` for the tests,
  which load its namespace explicitly before touching a result.

* The vendored core is refreshed to pick up `DfmNormalStochvol`, dfmtools'
  dynamic factor model with stochastic volatility in both error terms. Neither
  the sampler nor its header is compiled into this package -- `dfm_normal_stochvol.h`
  and `src/core/models/dfm_normal_stochvol.cpp` join `dfm_normal_gamma.h` and
  `dfm_normal_gamma.cpp` on `tools/update-bayests-core.R`'s skip list, for the
  same reason: nothing here includes either, since the dynamic factor models are
  dfmtools' rather than this package's. `DfmNormalStochvolInitial`,
  `DfmNormalStochvolInput`, `DfmNormalStochvolDraws` and
  `DfmNormalStochvolInput::validate()` are compiled in all the same, because they
  live in `inputs.h`, `results.h` and `inputs.cpp`, which every model shares and
  which are copied whole.

  Three shared files moved as part of that refresh, and none of them changes a
  draw of any model here. `chan_jeliazkov_2009`'s `precision_of` gained a
  diagonal fast path -- a scan and a reciprocal in place of a dense `inv_sympd`
  -- which every time varying parameter model in this package reaches with a
  diagonal state covariance and returns the same numbers from; `dfm_support.h`'s
  `initial_state_covariance` gained a second accepted shape that this package's
  models never pass; and `DfmNormalGammaInput::validate()`'s checks moved into a
  function two dynamic factor models now share, unchanged. Verified against the
  BayesTS source tree: the five touched files vendored here are byte-identical
  to upstream, `inst/COPYRIGHTS` needed no edit, and no vendored source collided
  with one of this package's own. Both additions to `model_support.h` and
  `chan_jeliazkov_2009` are new functions and a new accepted argument shape
  respectively, neither called by anything this package's models reach, so what
  every `VarNormalStochvol`, `VarTvpStochvol` and `VarTvpWishart` model already
  drew is untouched; a `VarTvpWishart` model drawn end to end after the refresh
  confirms the touched translation unit still compiles into a working sampler.
  This package has no test suite of its own; BayesTS's own golden fingerprint
  harness, the only regression check either side of the vendoring boundary has,
  passed 145 of 145 before this was vendored.

* Dynamic factor models have moved to [dfmtools](https://github.com/franzmohr/dfmtools). `create_dfmodel`, `dfmpost`, `dfm`, `add_priors.dfmodel`, `add_initial_values.dfmodel`, `plot.dfm`, `summary.dfm`, `thin.dfm` and the `bem_dfmdata` data set are gone from this package, as is the C++ helper `.post_lambda` behind `dfmpost`. That package now implements the same model through the vendored BayesTS core rather than in R, so its posterior simulation is faster and no longer limited to one factor and a transition of order at most two -- see its NEWS for what changed about the numbers. Nothing else here referenced them; the one internal mention was the list of classes `combine_models` accepts. The vendored core still carries the sampler's sources, `src/core/models/dfm_normal_gamma.cpp` among them, because `src/core/` is a verbatim mirror of upstream and is not curated per package -- it already carried `chan_jeliazkov_2009`, which nothing here calls either.
* Added `VecKlgs2010`, the vendored core's non-SUR implementation of the cointegration sampler of Koop, León-González and Strachan (2010), with the R entry points `.VecKlgs2010Coefficients`, `.VecKlgs2010Forecasts` and `.VecKlgs2010LogLik`. It draws the same posterior as `VecNormalWishart` -- from the same seed the two agree to 1.4e-14 on `a`, 2.4e-15 on `beta`, 3.2e-14 on the error precision and 6.0e-14 on the pointwise log likelihood -- and does it by exploiting the fact that a VEC's `k` equations share their regressors: the SUR design matrix is `kron(W_x, I_k)`, so the posterior precision of the coefficient block is `kron(W_x' W_x, Sigma^-1)` and can be formed from the compact regressors without building the SUR matrix at all. On a three-variable VEC of level order four over 160 periods and 1000 draws that is 0.05 s against 0.41 s. It reads `data$train$x`, which `create_bvecmodel` already assembles alongside `data$train$z`, so no model object has to change; what it does not offer is variable selection, since both schemes act on the columns of the matrix it declines to build, and it rejects either scheme with a message pointing at `VecNormalWishart`. Nothing in this package dispatches to it yet -- `draw_posterior.bvecmodel` still calls `VecNormalWishart` -- so no existing result changes.
* `add_priors` accepts a VEC model with time varying parameters, where it used to stop with "TVP priors need to be implemented." The cointegration vectors are then a state path rather than a draw from a cointegration space prior, so `coint` takes `rho`, the autocorrelation of that path, in place of the `v_i` and `p_tau_i` that describe the space's shrinkage and central location; `coef` takes `shape` and `rate` for the state error variances, as it already did for a TVP VAR. The prior on the cointegration state before the sample is the state equation's own stationary distribution, `N(0, I / (1 - rho^2))`, which is what makes it proper -- at `rho = 1` the path is a random walk whose variance grows without bound and which, `beta` being identified only up to scale, has nothing to pull it back. The loadings carry the compensating scale: only the product `alpha beta'` is identified, so their prior variance is shrunk by `1 - rho^2`, leaving the product on the scale `coef$v_i` asks for. Two things in the unreachable code behind that stop would not have worked as written and are fixed: the cointegration prior precision was named `v_i` where every reader of a normal prior in this package expects `v_inv`, and the prior precisions were built with `Matrix::Diagonal`, which the samplers cannot read -- they take a dense matrix, and an S4 sparse one fails the conversion with "Not a matrix." at the first draw. The prior precision of the non-loading coefficients now follows `coef$v_i` rather than being fixed at 1, as it already did for a TVP VAR and for a constant VEC. Nothing changes for a VEC with constant cointegration parameters.
* Added the R entry points for the five vendored VEC samplers that had none -- `VecNormalGamma`, `VecNormalStochvol`, `VecTvpWishart`, `VecTvpGamma` and `VecTvpStochvol` -- alongside the `VecNormalWishart` ones. Each is the same translation layer between the R model object and the core's structs that the six VAR models have, in `src/<Model>.cpp`, with `.<Model>Coefficients`, `.<Model>Forecasts` and `.<Model>LogLik`. What differs between them is which parts of the posterior each entry point needs: `beta` is read for all three, since without it `a` carries only the loadings and both the forecast and the log likelihood rebuild the loadings' regressors from the cointegration matrix; the time varying models slice `a` and `beta` to the last in-sample period to forecast from, by their own widths, neither of which is the column count of the forecast design; and the precision is read whole or at one period depending on whether it moves with time, which for the gamma models depends on whether there is a covariance block to move it. No draws change: nothing called these before, because they did not exist.
* The C++ VEC sampler `src/bvecalg.cpp`, the R one in `R/algo_bvectvp.R` and the `bvecpost` entry point are removed. No posterior simulation for VEC models is implemented in this package any more; all of it is the vendored core, as it already was for VAR models.
* **Draws of every time varying parameter model change**, and the previous ones were wrong by a period. `kalman_durbin_koopman_2002` returns `T+1` columns, of which the first `T` are the states the `T` observations load on and the last is the transition applied once past the end of the sample, informed by nothing. Every sampler kept the *last* `T` instead and paired them with the regressors of the first `T`, so the coefficient path was shifted one period ahead of the data. Within the same iteration that put the wrong period's coefficients into the residuals, which inflated the posterior scale of `Sigma`; into the state innovation variance, which saw a two-period jump as one increment, never saw the first increment at all, and counted a step drawn from the prior in its place, biasing it up; into the draw of the initial state, which used the wrong variance for the gap; and into the BVS likelihood ratio, which scored its candidates against the shifted path. On the way out, every reported coefficient path was a period ahead of the data, and its last period -- the one `add_posterior_forecasts` starts from -- carried a draw the data had never touched. This affects `VarTvpGamma`, `VarTvpWishart` and `VarTvpStochvol` in every configuration, and the five vendored VEC samplers that have no binding here. The size of the change scales with the state innovation variance, so a tightly shrunk path moves little and a freely moving one moves a lot. Constant coefficient models are untouched: they never call the smoother. The exported R function is unchanged -- it was always its callers that were wrong, and its documentation now says which column belongs to which period. Fixed in BayesTS; see its CHANGELOG for how the alignment was verified against a closed-form posterior.
* `add_posterior_forecasts` works for structural `VarTvpStochvol` models, where it used to fail with a message about the regressors and the posterior disagreeing on the number of coefficients. Three of the six vendored forecast routines never split the `K(K-1)/2` contemporaneous coefficients off the end of `a` and never applied `A_0^{-1}` to the simulated path; `VarTvpStochvol` is the only one of the three a structural model in this package can reach, since `create_bvarmodel` already refuses the other two error specifications. The forecast design this package builds carries no columns for those coefficients, so the mismatch was always an error and never a silently wrong path -- no forecast anyone obtained changes, there is simply now one where there was none. Fixed in BayesTS.
* The vendored core refuses a structural model whose error covariance is unrestricted -- a Wishart prior on the error precision, or a covariance block, `Psi` being a second contemporaneous matrix doing `A_0`'s job. `A_0` is unit lower triangular with `K(K-1)/2` free elements and the data determine only `A_0^{-1} \Sigma A_0^{-T}`, which has `K(K+1)/2`, so an unrestricted `\Sigma` leaves a `K(K-1)/2` dimensional set of parameters fitting identically and a draw of `A_0` carries no information. Nothing visible changes here: `create_bvarmodel` has always refused `structural = TRUE` with `error` in `wishart`, `gamma+covar` and `sv+covar`, and the two rules agree case for case. The check is a backstop for a model list assembled by hand or by another package linking the core.
* Added the vendored samplers for five more vector error correction models -- `VecNormalGamma`, `VecNormalStochvol`, `VecTvpGamma`, `VecTvpWishart` and `VecTvpStochvol` -- alongside `VecNormalWishart`. Like it, they are compiled into the shared object but have no binding and no R entry point, so this package's VEC models still run through `src/bvecalg.cpp` and `R/algo_bvectvp.R`. The three time varying ones are the port of `.bvectvpalg`, drawing the cointegration vectors as a state path with the same simulation smoother the coefficients use.
* Draws are otherwise unchanged for every model and configuration this package can reach. The refresh also collapsed the contemporaneous split, which three samplers had written out inline, into one pair of functions in the core, and moved the inversion of `A_0` out of the forecast horizon loop -- it does not depend on the horizon and was being redone `h` times per draw. Both were verified against a recorded fingerprint of every sampler before and after, byte for byte.

* The vendored core gained `chan_jeliazkov_2009`, the precision based alternative to `kalman_durbin_koopman_2002` for the same conditional posterior: the state path is one Gaussian vector whose precision is block tridiagonal, drawn in a single pass over a block banded Cholesky. It is compiled into the shared object but has no binding and no R entry point, and nothing in this package uses it -- the samplers still use the simulation smoother, because on the shapes a model here runs it is about twice as fast. The ratio is flat in both `T` and `M`: both algorithms are `O(T M^3)` and neither forms a `TM x TM` matrix, so there was no asymptotic advantage to be had, and per period the smoother does two matrix products where the precision based one does a Cholesky, a triangular solve and a symmetric product. What remains unexploited is the structure inside the blocks -- with a random walk the off-diagonal block is `-\Sigma_v^{-1}`, diagonal in `VarTvpWishart` and `VarTvpGamma`, and its triangular factor comes out exactly lower triangular, so a good half of the work is on known zeros. It is here as an independent implementation to validate the smoother against, and the same approach is already in use where it does win -- the stochastic volatility draw, where the state is scalar and the band is tridiagonal. Added in BayesTS; see its CHANGELOG for the measurements.
* `kalman_durbin_koopman_2002` no longer decomposes the same matrix once per period. Each of `sigma_u`, `sigma_v` and `B` may be supplied as one matrix that holds for every period or as a stack of one per period, and both forms remain supported in any combination -- a time varying error covariance, state variance and transition are all still available. What changed is that the constant form is no longer replicated into `T` copies of itself before the loops, which had left them unable to see that the blocks were identical, so a constant covariance was eigendecomposed `T` times. It is decomposed once. Draws are bit-identical, not merely unchanged to a rounding error: verified against the previous implementation over `T` from 2 to 120, `K` of 1 and 3, and all six combinations of constant and time varying arguments, and from R for the constant, stacked and time varying cases. 1.9 times faster at `T = 89`, `K = 3`, `M = 21` with constant covariances; the gain grows with `T` and vanishes when every argument really is time varying. This also speeds up `VarTvpGamma`, `VarTvpWishart` and `VarTvpStochvol`, which hand it constant covariances. Fixed in BayesTS.
* The exported R function `kalman_dk` is replaced by `kalman_durbin_koopman_2002`, which calls the vendored core instead of holding its own copy of the simulation smoother. `src/kalman_dk.cpp` was a third copy of the algorithm -- the same numerics as the core routine line for line, differing only in taking its arguments by value -- so a draw from a given seed is bit-identical to what `kalman_dk` returned, and no posterior of any model changes. There is no speed difference either; this removes the duplicate, nothing else. The example in the documentation was rewritten, since it called `gen_var` and reached for `temp$data$Y` and `temp$data$SUR`.
* The C++ entry point `bvartools::kalman_dk` no longer exists. It was added in 0.2.2 and is a removal rather than a tidy-up: a package that called it has to link `src/core/algorithms/kalman_durbin_koopman_2002.cpp` instead. The generated interface wraps every call in an `Rcpp::RNGScope`, which rewinds the random number stream for a caller that holds the RNG across calls -- the defect this release documents three instances of -- so the replacement deliberately does not generate one.
* The exported R functions `stochvol_ksc1998` and `stochvol_ocsn2007` are replaced by `stochvol_ksc_1998` and `stochvol_ocsn_2007`, which call the vendored core instead of reimplementing it. Between 8 and 93 times faster over `T` from 100 to 2000 at `K = 3` -- the R versions materialised the `T x T` posterior precision of the log-volatility and factorised it, where the core exploits the fact that it is tridiagonal and factorises it in `O(T)`. The two implementations of each algorithm that had to be kept in step by hand are now one, which is how `stochvol_ksc1998` came to be the only one of the four missing the mixture-underflow and input-validation fixes of the previous release. Draws from a given seed differ, since the compiled draw consumes the random number stream in a different order; the posteriors agree -- over 300 sweeps at `T = 300`, `K = 3` the root mean squared error against a known log-volatility path is 0.424 against 0.427 for Kim, Shephard and Chib and 0.426 against 0.434 for Omori, Chib, Shephard and Nakajima, and the two posterior means differ by at most 0.27 where the path itself has a standard deviation of 1.87. `set.seed` still reaches the draw. Bad arguments are now reported by the core and name themselves, as in `stochvol_ksc_1998: 'sigma' must have 3 elements, got 2`.
* The stochastic volatility draw of the vendored core now factorises the posterior precision of the log-volatility as the banded matrix it is. It is tridiagonal -- the random walk contributes the two bands of `D'D` and the normal mixture only a diagonal -- and the Cholesky factor of a symmetric tridiagonal is bidiagonal, so one sweep over the periods replaces a dense factorisation that cost `O(T^3)` in time and `O(T^2)` in memory. Draws are unchanged bar rounding (relative difference 1.8e-15 from the same seed). This also speeds up `VarNormalStochvol`, `VarTvpStochvol` and the VEC sampler with stochastic volatility, which reach the same routine.
* Added `stochvol_ksc_1998`, the seven-component mixture of Kim, Shephard and Chib (1998), as a counterpart to the ten-component `stochvol_ocsn_2007`. No sampler in this package uses it; it is available for one written in R.
* `.bvectvpalg` calls `stochvol_ksc_1998` in place of `stochvol_ksc1998`, so its draws from a given seed change. That sampler cannot currently be reached from R -- `add_priors` refuses a VEC model with time varying parameters -- so nothing observable changes.
* Migrated `bvarpost` into `add_posterior_coefficients.bvarmodel`
* Added the core layer of BayesTS -- the samplers' declarations in `inst/include/bayests` and their numerics in `src/core` -- as vendored source. See `src/core/VENDORED.md`.
* All six `bvarmodel` algorithms -- `VarNormalWishart`, `VarNormalGamma`, `VarNormalStochvol`, `VarTvpGamma`, `VarTvpWishart` and `VarTvpStochvol` -- are now the samplers from the vendored core layer, and each model's three entry points are consolidated into one source file, `src/<Model>.cpp`. No posterior simulation is implemented in this package any more. Draws agree with the previous implementations to a rounding error (worst relative difference 5.6e-14 across eighteen model configurations; all inclusion indicators identical).
* BVS selection of the contemporaneous error coefficients was not selecting in `VarTvpGamma` and `VarTvpStochvol`. Those coefficients are a path, one column per period, but the candidate the likelihood ratio scored zeroed the coefficient in the first period alone rather than across the sample -- a single linear index into the matrix instead of a whole row. The two candidates being compared therefore differed in one period out of `T`, the ratio between them was correspondingly close to zero, and inclusion won essentially every time: at `T = 24` all three free coefficients stayed in for all of a recorded 80 draws, against 26 of 240 once the candidate spans the sample. The mask applied to the draws that are reported was never affected, so the coefficients themselves were consistent with the indicators -- it was the indicators that carried no information. Fixed in BayesTS. Draws of `VarTvpGamma` and `VarTvpStochvol` with both a covariance block and variable selection change; no other model or configuration is affected.
* The two variable selection schemes of the vendored core now live in `src/core/algorithms/ssvs.cpp` and `src/core/algorithms/bvs.cpp`, one file per scheme, in place of thirteen near-identical copies spread across the six samplers. Draws from a given seed are unchanged for every model and configuration.
* `VarTvpWishart` draws change beyond a rounding error. Its sampler built the block diagonal error precision once, before the first iteration, and then kept using it although the precision is redrawn every iteration, so the coefficient smoother conditioned on the initial error covariance for the whole chain. Fixed in BayesTS.
* The `VarNormalWishart` algorithm for objects of class `bvarmodel` is now the sampler from the vendored core layer. `bvarpost`, `add_posterior_forecasts` and `add_posterior_loglik` are unchanged and still dispatch to `.VarNormalWishartCoefficients`, `.VarNormalWishartForecasts` and `.VarNormalWishartLogLik`, which have become bindings over `bayests::VarNormalWishartSampler`; the R model object and the posterior it returns keep their shape. Draws agree with the previous implementation to a rounding error (relative differences below 1e-13, inclusion indicators identical), which comes from the core reusing one Cholesky factor for the posterior mean and the draw and inverting the Wishart scale with `inv_sympd`. Inconsistent input is now reported by name -- "initial error precision must be 3x3, got 2x2" -- rather than as a matrix dimension mismatch.
* The simulation smoother of Durbin and Koopman (2002) now comes from that core layer, and the TVP samplers call it directly instead of through the package's own C++ interface. The interface wrapped every call in an `Rcpp::RNGScope`, which rewound the random number stream and handed the same numbers out twice per iteration, so posterior draws of TVP models change from the second iteration on. The C++ entry point `bvartools::kalman_durbin_koopman_2002` no longer exists; it was never available from R.
* Fixed the same random number stream rewind in `post_bvs`, which called `bvartools::sur_const_to_tvp` for time varying parameters, and in the VEC sampler, which called `bvartools::stochvol_ocsn2007_internal` for stochastic volatility. Both now call the function directly. Draws of the affected models change. The `bvartools::sur_const_to_tvp` wrapper itself is unchanged and remains available to other packages; `bvartools::stochvol_ocsn2007_internal` has since been removed, see the next entry.
* The stochastic volatility draw of Omori, Chib, Shephard and Nakajima (2007) now also comes from the vendored core layer, as `stochvol_ocsn_2007()`, and replaces `src/stochvol_ocsn2007_internal.cpp`. It validates its arguments instead of indexing them on trust, and draws the mixture indicator in logs, so an observation far out in the tails of all ten mixture components no longer normalises to `NaN` and selects a component that does not exist. The C++ entry point `bvartools::stochvol_ocsn2007_internal` no longer exists; it was never available from R. Draws of the VEC sampler with stochastic volatility change by a rounding error, because the posterior mean is now formed from the same Cholesky factor as the draw rather than by a separate LU solve.
* The exported R function `stochvol_ocsn2007` is a separate implementation and received the same three fixes. It validates `y`, `h`, `sigma`, `h_init` and `constant` rather than reporting a bad argument as a failed matrix factorisation or, for `h_init`, not at all; it draws the mixture indicator in logs, which stops an observation far out in the tails of all ten components from returning a matrix of `NA`; and it reuses one Cholesky factor for the posterior mean and the draw. Draws from a given seed are unchanged up to a rounding error. The documentation of the mixture, previously described as having seven components, now says ten.
* Fixed missing transformations for structural modesl in `irf.bvar` and `fevd.bvar`.
* Added function `choose_best_model`.
* Added function `create_first_difference_matrix`.
* Added function `create_second_difference_matrix`.
* Added functions `create_bvarmodel`, `create_bvecmodel` and `create_dfmodel` to replace `gen_var`, `gen_vec` and `gen_dfm`, respectively, in the future.
* Added `generate_lower_block_diagonal` for faster simulation of autocorrelated TVP coefficients.
* Updated functions for dynamic factor models to the revised design pattern.
* Added a warning message that the functionality for dynamic factor models will be exported to a separate package with package in the future.
* Add a generic function `mean_absolute_forecast_error` for MAFE calculation.
* Add a generic function `prediction_matrix` which helps with forecast generation.
* Add a generic function `selection_criteria` to calculate information criteria for model selection.
* Add a generic function `forecast_errors` to generate forecast errors.
* Add a generic function `add_initial_values` to separate initial value generation from prior specification.
* Added `gen_artificial_vec` to generate artificial data sets for algorithm and model testing.
* Added `gen_artificial_var` to generate artificial data sets for algorithm and model testing.
* Added `coint_kls2010_reparameterise_two` for more convenient data transformation.
* Added `coint_prepare_sur_data` for more convenient input data preparation for cointegration simulation.
* `plot.bvar` and `plot.bvec` allow to specify whether a horizontal line should be added or not.
* Added convenience function `covar_vector_to_matrix`.
* Added convenience function `sur_const_to_tvp`.
* General updates in documentations including update of `Rcpp` dependency in DESCRIPTION file to version 1.0.12.
* Added `post_gamma_state_variance` for posterior simulation of constant error variances of the state equation.
* Added `post_gamma_measurement_variance` for posterior simulation of constant error variances of the measurement equation.
* Renamed `.prep_covar_data` to `covar_prepare_data` and made it visible in R and also callable from C++.

# bvartools 0.2.4

* Using an updated version of `Rcpp` to address an issue with `Rcpp::stop`.
* `stochvol_ocsn2007` can handle multi-column input.
* `stochvol_ksc1998` can handle multi-column input.
* Added `post_normal_covar_tvp` for posterior simulation of time varying, lower triangular covariance matrices.
* Added `post_normal_covar_const` for posterior simulation of constant, lower triangular covariance matrices.

# bvartools 0.2.3

* Fixed alias issue resulting from use of `roxygen2`.
* Made `kalman_dk` callable from C++.
* Stochastic volatility algorithms allow to set the offsetting constant manually.
* Changed `stoch_vol` to a wrapper for `stochvol_ksc1998`.
* Added stochastic volatility algorithm of Kim et al. (1998) in a separate function `stochvol_ksc1998`.
* Added stochastic volatility algorithm of Omori et al. (2007) in function `stochvol_ocsn2007`.
* Fixed bug with detection of deterministic terms in `bvar`.
* Implemented recursive iterations for forecasts in C++.
* Replaced erroneous `|` in C++ sampling functions by `||`.

# bvartools 0.2.2

* Addressed CRAN NOTE on CITATION file
* Addressed the CRAN NOTE "Specified C++11: please drop specification unless essential" by dropping the specification from "src/Makevars"
* Improved the treatment of `bvar` and `bvec` objects if Gibbs sampler fails.
* Fix erroneous SUR-matrix generation for VEC models with r = 0 in `.bvecalg`.
* Fix bug in `.bvecalg` and `.bvectvpalg` with the storing of posterior draws of beta.
* Fix bug of `predict.bvar`, which could not handle only VARX models with contemporaneous exogenous variables only.
* Model plot functions support boxplots.
* Fix typos in documentation.

# bvartools 0.2.1

* Added functionality for the simulation of models with time varying parameters, both for VAR and VEC models.
* Added functionality for the simulation of models with stochastic volatility, both for VAR and VEC models.
* Added a plot function for classes `bvar` and `bvec` for visual inspection of posterior draws.
* Changed the generation of the output object in the Gibbs sampler functions `bvaralg` and `bvecalg` to make them more stable for especially large output.
* Changed `draw_posterior` to a generic function and added the corresponding methods for BVAR, BVEC and DFM input.
* Changed `irf` and `fevd` to generic functions.
* Corrected typos in documentation.
* `thin_posterior` methods were renamed to `thin` and are now methods of `coda::thin`.
* Function `irf` allows to specify the size of a shock.
* Fixed a bug in `ssvs_prior` concerning BVEC models.
* Fixed a bug with the prior in the BVEC algorithm.

# bvartools 0.2.0

* Changed `thin_posterior` to a generic function and added methods for BVAR, BVEC and dynamic factor model input.
* Changed `add_prior` to a generic function and added methods for BVAR, BVEC and dynamic factor model input.
* Added funcionality to estimate dynamic factor models (DFM).
* `predict` requires to specify an object of class `ts` as input for argument `exogen`.
* Additioal argument checks for `add_priors` methods.
* Updated documentation in `minnesota_prior` and for `add_prior` methods.
* Using \doi instead of \url in documentation

# bvartools 0.1.0

* Omitted package `Matrix` from "Imports"" in DESCRIPTION, which caused a note in version 0.0.3.
* Added function `bvarpost` for posterior simulation of BVAR models.
* Added function `bvecpost` for posterior simulation of BVEC models.
* Added function `draw_posterior` for estimation of multiple models.
* Fixed erroneous calculation of structural forecast error variance decompositions.
* More specification checks and increased robustness against erroneous model specificaions.
* Function `fevd` calculates FEVDs based on means of posterior draws of FEVDs and not based on the means of the coefficient draws.
* Function `bvar` and `summary.bvar` can deal with inclusion parameters.
* Added funtion `add_priors` for easier construction of prior matrices for multiple models.
* `gen_var` and `gen_vec` can produce multiple models.
* Changed all argument names of `predict.bvar` to lower cases.

# bvartools 0.0.3

* Changed all argument names of `post_normal`, `post_normal_sur`, `post_coint_kls` and `post_coint_kls_sur` to lower case letters.
* Replaced output element in function `ssvs` from `V_i` to `v_i`.
* Refined function `minnesota_prior` and added additional functionaliy.
* Fixed error message when creating seasonal dummies with `gen_var` and `gen_vec`.
* New data set `us_macrodata`.
* Added additional checks in `gen_vec`.
* Added functions `inclusion_prior` for the calculation of inclusion probability priors as used by `bvs` and `ssvs`.
* Added `summary` functions.
* Fixed conversion and collection of exogenous regressors in `bvec_to_bvar`.
* Fixed detection of deterministic terms in `bvec_to_bvar`.
* Updated documentation in `kalman_dk`.
* `irf` contains a new argument `keep_draws`.
* Additional checks in `post_normal`, `post_normal_sur`, `post_coint_kls` and `post_coint_kls_sur`.
* Adapt vignette `bvec`.
* Added `loglik_normal` for the calculation of a multivariate normal log-likelihood.

# bvartools 0.0.2

* Updated vignette `ssvs` after the introduction of function `ssvs_prior`.
* Added `ssvs_prior` for the calculation of prior matrices for the SSVS algorithm.
* Added `minnesota_prior` for the calculation of the Minnesota prior.
* Use unsigned integers for indices in Cpp code to address warnings during installation.
* Better error handling in `irf`.
* In `post_coint_kls_sur` the prior matrix `g_i` can be time varying.
* `bvar` and `predict` also work only with deterministic terms, i.e. p can be zero.
* Use SVD to obtain a draw of beta in `post_coint_kls` and `post_coint_kls_sur`.
* `predict` allows for p = 1.
* Add legend to `plot.bvarfevd`.

# bvartools 0.0.1

* Initial release
