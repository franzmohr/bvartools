# bvartools 1.0.0

* **The forecasts of a discounted model are `coda::mcmc` draws, and an
  estimated discounted model can be written to HDF5.** `add_posterior_forecasts()`
  returned the i.i.d. forecast draws of `VarTvpDiscount` and `VecTvpDiscount` as
  a plain matrix, without the start, end and thinning interval every step after
  it reads, so `add_forecast_errors()` stopped with "non-numeric argument to
  mathematical function" and `write_to_hdf5()` with "If a sample robj is not
  provided, both dtype and space have to be given". They are now labelled as
  draws 1 to `iterations`, unthinned, which is what BayesTS writes for them.
  `write_to_hdf5()` also stopped on the per-period posterior itself, and would
  have left out `posterior$u_sigma` and `posterior$df`; it now writes the
  closed-form blocks without the chain attributes, as BayesTS does, and
  `read_model_from_hdf5()` reads them back as the plain matrices they are.
  `selection_criteria()` on an expanding window of discounted models reports
  the last window's `LML`, which it dropped. The draws are unchanged: only
  their labels are new.

* **`plot()` draws the criteria of an expanding window.** `selection_criteria()`
  on an `expandingwindow` or a single model returns a `selcrit` that keeps the
  model's classes, and with no `plot.selcrit()` the call reached
  `plot.bvarmodel()` and stopped with "dim(X) must have positive length". The
  new method plots it as a `selcritlist` of one. `plot()` of either also
  accepts `criterion = "LML"`.

* **`create_external_forecast()` refuses forecasts at another frequency than
  the data.** Periods were rounded to the frequency of the models, so annual
  forecasts for 2020 and 2021 given to a quarterly model became forecasts of
  2020Q1 and 2021Q1, and an annual growth rate was scored as a quarterly one.
  The frequency is read off the spacing of the periods within a publication:
  periods further apart than one period of the data, or several that fall into
  the same period, now stop with a message saying which frequency they appear
  to be. Where every publication forecasts a single period, periods that all
  sit at the same position within the year are refused as well.

* **`VecNormalGamma` and `VecNormalStochvol` draw from the posterior of the
  prior they state, and their draws change.** Both use the cointegration space
  prior of Koop, Leon-Gonzalez and Strachan (2010), which scales the loadings'
  prior by G. For `VecNormalGamma`, G is the error covariance, so the prior is
  also a factor in the error precisions' posterior, and the vendored core left
  it out. `VecNormalWishart` has always included it. For `VecNormalStochvol`,
  G was re-averaged from the current volatilities after every draw, which made
  the prior a function of them that their own draw never saw. It is now fixed
  for the run, from the starting volatilities. *Draws change* for both models
  with `rank > 0`, in every configuration. With one variable,
  `VecNormalGamma` and `VecNormalWishart` given matching priors are the same
  model, and upstream's test finds them within 0.07% of each other on the
  posterior mean of the precision, against 2.6% before. `VecNormalWishart`
  and `VecKlgs2010` draws are unchanged.

* **`coint$g_i` sets the G of `VecNormalStochvol`.** `add_priors()` and
  `cointspace_prior()` take the inverse of G, the matrix the loadings' prior is
  scaled by, as diagonal elements, a full matrix or `"ml"` for the inverse of
  Johansen's error covariance, and store it as `priors$beta$g_inv`, the name
  a BayesTS model file uses. Left out, G still comes from the starting
  log-volatilities, so the prior depends on `add_initial_values()`; given, it
  does not. It is refused for every other model, which scales the loadings'
  prior by its error covariance.

* **The VEC models refuse a cointegration prior the sampler cannot honour.**
  The prior on the loadings has to be centred at zero and independent of the
  other coefficients, so a non-zero prior mean on the first `k * rank`
  coefficients, or a prior precision coupling them to the rest, now stops with
  a message naming the position. So do a negative `v_i` and a `p_tau_inv` that
  is not positive definite while `v_i` is positive. `add_priors()` and
  `cointspace_prior()` produce none of these, so only a prior edited by hand
  reaches them.

* **`bayests_files()` renames the exit status into place once it is written.**
  The script it runs wrote BayesTS's exit status straight into the file the
  session polls for, and a redirection creates that file before anything is in
  it, so a poll could find it empty and report a run that had succeeded as
  failed with "exit status NA". It happened once, on 21 September 2026, to two
  directories of a pass of `bayests loglik` that mostly skipped files that
  already had their log-likelihood, which finish well inside the polling
  interval. It did not happen again in 20 repetitions of that pass or in 180
  runs polled every 10 ms, so it is rare; the status now goes to a file of its
  own and is renamed into place, and a status that still cannot be read is read
  again for a quarter of a second before it counts, in case something else held
  the file open for a moment. A run that fails is still reported, with the
  status it exited with.

* **SSVS refuses three priors it could not honour.** The vendored core now
  checks, at every coefficient SSVS selects over, that the prior mean is zero,
  that the prior precision couples it to nothing else, and that `tau0` is
  smaller than `tau1`. Each of those files used to run: the coefficients were
  drawn under the prior as given while the inclusion indicators were scored as
  if both mixture components sat at zero, so the chain targeted neither model
  (George, Sun and Ni 2008, eq. 12). `add_priors()` produces a diagonal
  precision and `tau0 < tau1` from its defaults, so what reaches this is a
  non-zero mean: **`coef$const` together with SSVS and
  `varsel$exclude_det = FALSE`, the default, now stops** with a message naming
  the position, as does `coef$coint_var = TRUE` with SSVS over the first own
  lags it puts a mean of one on. Setting `varsel$exclude_det = TRUE` keeps the intercept out of the
  selection and runs. *Draws are unchanged* for every model still accepted.

* **`bvs` warns when its coefficient prior is too flat to select against.** An
  excluded coefficient is drawn from its prior and then scored against the
  data, so the flatter the prior the harder it is for anything to get back in,
  and inclusion probabilities pinned near zero describe the prior rather than
  the data (Korobilis 2013, section 3.1). The core reports a selected position
  whose conditional prior variance is 100 or more, and
  `add_posterior_coefficients()` raises it as an R warning for the seven
  constant-coefficient models that offer `bvs`. It was collected rather than
  raised where it arose: `Rf_warning()` may longjmp out of a running sampler.
  *Draws are unchanged.*

* The vendored BayesTS core is refreshed to upstream `cea124b`, which brings
  the two items above and the one below. *Draws are unchanged* for every model
  that does not use the new prior: upstream's fingerprint recording is
  identical for every file without it.

* **A test for time variation in TVP-VARs with stochastic volatility.**
  `add_priors()` takes `coef$omega_v` and `sigma$omega_v` in place of
  `shape`/`rate` for a model with `tvp = TRUE` and `error = "sv"` or
  `"sv+covar"`. Each puts a normal prior `N(0, omega_v)` on the signed standard
  deviation of a random walk's innovations -- the coefficients and the
  covariance coefficients for `coef`, the log-volatilities for `sigma` -- which
  is the non-centred parameterisation of Frühwirth-Schnatter and Wagner (2010).
  A coefficient or volatility that does not move is then a point inside the
  prior, so the Bayes factor for time variation is a Savage-Dickey density ratio
  that one run estimates (Chan 2018). Each block of the posterior drawn this way
  holds, beside `sigma` (still the state variance, now `omega^2`), the draws of
  `omega` and the log densities at zero `omega_log_zero`, per state, and
  `omega_log_zero_joint`, for the block; `?add_priors.bvarmodel` gives the
  formula. The blocks choose their prior one by one, `omega_v` is refused
  beside `shape`/`rate` and on every other model, and `write_to_hdf5()` and
  `read_model_from_hdf5()` carry the prior and the draws. `VarTvpStochvol` now
  also returns the core's warnings, which `add_posterior_coefficients()`
  raises.

* **`time_variation_test()` turns those draws into Bayes factors.** For a
  `bvarmodel` estimated under `omega_v` it reports, for every coefficient,
  covariance coefficient and log-volatility, the log Bayes factor in favour of
  time variation against a constant state, and one for each block jointly,
  each with a numerical standard error from batch means. The joint Bayes factor
  compares "every state of the block moves" with "none does", not "at least one
  moves", and `?time_variation_test.bvarmodel` says why the two can disagree.

* **The test covers TVP-VARs with a gamma error term too.** `coef$omega_v` is
  accepted for `tvp = TRUE` with `error = "gamma"` or `"gamma+covar"` as well,
  for the coefficients and the covariance coefficients, whose draws and
  `time_variation_test()` rows are those of the stochastic volatility model.
  The error precision of these models does not move, so `sigma$omega_v`
  remains for stochastic volatility, and is now refused on every other error
  term rather than accepted as a name and left unread.

* **And TVP-VEC models with stochastic volatility.** `add_priors()` on a
  `bvecmodel` with `tvp = TRUE` and `error = "sv"` or `"sv+covar"` takes
  `coef$omega_v` (the loadings, the other coefficients and the covariance
  coefficients) and `sigma$omega_v` (the log-volatilities) in place of
  `shape`/`rate`, and `time_variation_test()` has a `bvecmodel` method, which
  labels the loadings by the error correction term they load on. The
  cointegration space keeps its state equation, whose variance is fixed to pin
  down the scale of beta, so it has no prior to test against. The same
  `coef$omega_v` works for TVP-VEC models with `error = "gamma"` or
  `"gamma+covar"`, whose coefficients and covariance coefficients are tested
  as those of the stochastic volatility model are.

* **`transform_variables()` returns a vector series for a vector series**, where
  it returned a one-column matrix, and a single series takes a named or an
  unnamed code alike: a named code on a series without a column name used to
  stop with "subscript out of bounds". Its documentation said code 7 loses one
  leading observation; it loses two, being the difference of a growth rate, and
  always did.

* **Every plot method returns its input invisibly**, which `plot.bvarfevd()`,
  `plot.bvarirf()`, `plot.bvarprd()`, `plot.modellist()` and
  `plot_forecast_errors_by_period()` did not: they returned whatever their last
  drawing call did. **`write_to_hdf5()` on a 'modellist' or an
  'expandingwindow' returns the paths it wrote**, invisibly, as it already did
  for a single model.

* Every exported function documents its value. The help of `post_coint_kls()`
  and `post_coint_kls_sur()` said `Gamma` is a K x N matrix; it is a column
  vector, `vec(Gamma)` for the first, which is how `bvec()`'s example and the
  VEC vignette have always used it. `selection_criteria()` says what the median
  and the band of `RSFE` are -- the ones of `AFE` up to the interpolation
  between draws -- and the agent skill in `inst/agents` covers the discounted
  models and `LML`.

* `Depends` requires `R (>= 4.0.0)` rather than `(>= 3.5)`. `src/Makevars` has
  set `CXX_STD = CXX17` since the C++ core arrived, and 3.5 predates R's own
  requirement of a C++11 compiler, let alone a toolchain that honours a request
  for C++17 -- so the old floor was a claim nothing tested, the check matrix
  reaching back only to `oldrel-1`. dfmtools and fincond, which carry the same
  mismatch or depend on this one, move with it.

* **The two discounted models, `VarTvpDiscount` and `VecTvpDiscount`, can be set
  up, estimated, forecast, scored, written and read back.** They are reached
  with `algorithm = "discount"` of `create_bvarmodel()` and
  `create_bvecmodel()`, and `add_posterior_coefficients()` estimates them here
  like any other algorithm: the filter is part of the vendored BayesTS core, as
  every sampler in this package is. Estimating one consumes no random numbers,
  so a model estimated in R and the same model estimated by the `bayests`
  command line over a written file agree to the bit rather than merely in
  distribution. No sampler here changes, and no model that does not ask for the
  new algorithm is written differently.

    What comes back is a posterior rather than a chain, and the steps that
    follow know it: `posterior$a$mean`, `posterior$a$scale`, `posterior$a$cov`,
    `posterior$u_sigma$scale` and `posterior$df` hold one row per period and
    carry no `mcpar`, and there is no `coeffs` anywhere, because joining one
    draw per period would look like a sampled path and is not one.

    They are not samplers. The posterior is the matrix normal dynamic linear
    model of West & Harrison (1997, ch. 16) with the discounted Wishart of Uhlig
    (1997), closed form in one pass over the sample, so `burnin` must be 0 and
    `thin` 1 and `iterations` says only how many i.i.d. draws a forecast takes.
    Two discounts govern it, `delta_beta` for the coefficients and
    `delta_sigma` for the error covariance, both in `(0, 1]` and both a model in
    their own right at 1, where the quantity they govern does not move. A vector
    in either gives one model per value, as a vector in `p`, `s` or `r` does.

    What they buy, besides the speed, is that the sum of `/posterior/loglik` is
    the **exact** log marginal likelihood of the sample rather than an estimate
    of it. `selection_criteria()` reports it as `LML`, and it is the one
    criterion they carry: there is no chain to estimate an effective number of
    parameters from, and the marginal likelihood has already paid for the
    complexity a count of parameters would charge for. `choose_best_model()`
    maximises it, so a grid over the rank, the lag order, the cointegration
    matrix or the two discounts is compared without a chain being run for any of
    it.

    Three things about the file differ from every other model here, and each is
    a different model rather than a spelling. The coefficient prior is a matrix
    normal at `/priors/a/mean` and `/priors/a/cov` -- an n_design x k mean and
    the n_design x n_design regressor side of a covariance -- and not the
    `/priors/a/mu` and `/priors/a/v_inv` of a sampler, which a discounted model
    reads as no prior at all. `add_priors()` therefore takes `coef$v_i`,
    `coef$v_i_det`, `coef$v_i_alpha` and `coef$const`, refuses `coef$shape` and
    `coef$rate`, which are the prior of state variances this model does not
    have, and refuses `coef$v_i = 0`, a precision with no covariance to write.
    The SUR matrix `/data/train/z` is not written, the filter running against
    the compact `/data/train/x` instead. And a discounted VEC **conditions on a
    fixed cointegration matrix** rather than drawing one:
    `add_initial_values()` puts Johansen's estimate at `/initial/beta`, or the
    space given in its new `beta` argument, and there is no cointegration space
    prior to specify. What drifts is the adjustment to the long-run relations
    and not the relations themselves, which is a different question from the one
    `algorithm = "KLGS2010"` answers rather than a cheaper way of answering the
    same one.

    `open_models()` reports `delta_beta` and `delta_sigma` in the manifest, at 1
    for every model that has no discounts, which is what 1 means.

* **`read_model_from_hdf5()` reads a posterior block of one row as one row.**
  hdf5r drops a dimension of size one, and `as.matrix()` then made a column of
  what the file holds as a row, so such a block came back transposed: a single
  draw over many columns read as many draws of one column. Nothing a sampler
  writes has one draw, which is why this went unseen -- the discounted models'
  `/posterior/loglik` and their one-column `/posterior/beta/coeffs` are exactly
  that. A block with no `mcpar` attributes is now also read as the plain matrix
  it is rather than labelled as a chain, because a closed form's columns are
  periods and its rows are not draws.

* **Vendored BayesTS core refreshed to BayesTS 0.3.0 plus one commit,
  `6fe91d2`.** **Draws are unchanged** for every
  sampler here, and nothing this package compiles behaves differently: the four
  headers that changed -- `bayests/inputs.h`, `priors.h`, `results.h` and
  `spec.h` -- gain declarations and nothing else, and no vendored source reads
  what they declare. The fix itself does not reach this package. It is to
  `core/models/factor_score.h`, the filter a *factor* model's forecast is scored
  by, which is one of the sources `src/core/VENDORED.md` records as not copied;
  a VAR or a VEC is scored by its own pointwise log likelihood over the scored
  periods, which needs no filter. dfmtools, which does vendor that file, carries
  the fix.

    What 0.3.0 adds is `VarTvpDiscount` and `VecTvpDiscount`, two models with a
    closed-form posterior rather than a chain, and **both are now vendored**:
    `core/models/var_tvp_discount.cpp`, `core/models/vec_tvp_discount.cpp` and
    the `core/models/discount_support.h` they share. An earlier refresh held
    them back because nothing here could reach them; `src/VarTvpDiscount.cpp`
    and `src/VecTvpDiscount.cpp` now can, so their entries are out of the
    refresh script's `skip` and `src/core/VENDORED.md` describes them under a
    section of their own rather than under *Not copied*.

    The one commit past the release is a fix to `VarTvpDiscount`, found by
    wiring its score up here: `predictive_log_density()` read the rows of
    `/data/forecast/x` as they arrived, and the lagged endogenous blocks of a
    horizon past the first hold a placeholder, which a forecast overwrites as it
    simulates and that recursion does not. The first horizon was scored
    correctly and every one after it came back as `NaN`. It now fills those
    blocks from the realised values, as the eight sampling VARs beside it always
    have. `VecTvpDiscount` was never affected.

    `src/core/VENDORED.md` now also records which upstream commit the copy is
    at, which it never has. Nothing but a reader enforces that paragraph -- the
    refresh script compares files and `inst/COPYRIGHTS`, not prose -- so it says
    to check it against the upstream log during a refresh.

* **`selection_criteria()` has a default method.** A model of a class this
  package has no method for -- a dynamic factor model of dfmtools, say -- now
  gets `LL`, `WAIC`, `LOOIC` from `posterior$loglik` and `LPL` from
  `posterior$forecast$loglik`, which is everything the draws alone support.
  `choose_best_model()` ranks the result and `print()` shows it, both unchanged.

  What the default cannot report is `AIC`, `BIC` and `HQ`, which charge a model
  for its size and so need a count of its free parameters, and `FE`, `AFE` and
  `RSFE`, which need the variables named and paired up. Those are properties of
  a model rather than of its draws, and a class that has them should write its
  own method. The periods of the `LPL` terms are numbered from one for the same
  reason: this method does not know where a model keeps its sample, so it cannot
  say when the scored periods were.

  The point is that one implementation serves every class. WAIC and LOOIC in
  particular need nothing but a pointwise log-likelihood, and a package that has
  one should not have to write them again to report them.

* **Dynamic multipliers.** `multipliers()` returns the response of an endogenous
  variable to a change in a weakly exogenous one, with methods for a
  `bvarmodel` and a `bvecmodel`. The change is held from period zero on by
  default, or confined to period zero with `type = "transitory"`. An error
  correction model is put into its levels form with `vec_to_var()` first, so
  the multipliers are statements about the levels and the long-run relations
  enter them. What comes back has the class an impulse response has, so it
  plots the same way. This is the quantity the country models of a global VAR
  are read for, where the foreign block is weakly exogenous.

* **`selection_criteria()` reports the score of a forecast as `LPL`.** A model
  that carries `posterior$forecast$loglik` -- the log predictive density of each
  period its horizon realised, written by BayesTS against `data$test$y` -- now
  gets the criterion `"LPL"`, the log of the mean of the draws of each period's
  density, summed over the periods. `choose_best_model()` takes it as before and
  `print()` shows it beside the in-sample criteria.

  It is the criterion an expanding window exercise already reported, from the
  densities of `add_predictive_loglik()`, and the two now go through one
  function: the same densities give the same entry, band and numerical standard
  error alike, so a comparison of a scored model with a scored window is a
  comparison of the same quantity. A single model's densities are those of the
  horizons of one forecast, each conditioning on the periods realised before it;
  an expanding window's are one per window. Both are one step ahead and both sum
  to the log predictive likelihood of the stretch they cover.

* **A model carries what it was scored against, in `data$test$y`.**
  `add_forecast_errors()` puts the periods of the horizon it took the errors
  against into the model, and `write_to_hdf5()` writes them to `/data/test/y`,
  beside `/data/train/y`. `test_sample` now defaults to `NULL`, in which case
  those values are used, so a model read back from a file is scored without the
  sample being supplied a second time.

  This is what lets an expanding window exercise be scored a file at a time. A
  window used to need its test data passed in from outside, which for a folder
  of models means holding the series in the session and knowing which periods
  belong to which window; now each file carries the periods it is to be judged
  by. `vec_to_var()` takes them across with the forecast data, the levels being
  what a VEC model is scored in either way.

  BayesTS 0.3.0 reads the same dataset into every model's input, so `bayests
  check` no longer reports it as one the model never reads, and the predictive
  log-likelihood that will fill `posterior$forecast$loglik` has its observations
  waiting for it. Nothing in BayesTS computes from it yet.

* **The forecasts move to `posterior$forecast$forecasts`, with the forecast
  errors beside them at `posterior$forecast$errors`.** `posterior$forecast` is
  a group now rather than a matrix of draws, and `posterior$forecast_errors` is
  gone. **This breaks code that reads either of the old names**, in a model
  object as much as in a file written with `write_to_hdf5()`; the accessors
  `predict()`, `get_forecast_errors()` and `selection_criteria()` are unchanged
  and are the way to reach both without naming a path.

  The group is the place for everything the forecast periods produce, which is
  about to include the log predictive density of what those periods realised.
  One matrix could not hold three things, and the errors were already spelling
  the group with an underscore. The members are named after what they hold
  rather than one of them being `draws`, since all of them are draws.
  `posterior$loglik`, the in-sample pointwise log-likelihood, does not move: it
  evaluates each observation of the sample under states that have already seen
  it, which is a different statistic from a forecast score rather than the same
  one over other periods.

  It is the layout of the model file as well -- `/posterior/forecast/forecasts`
  and `/posterior/forecast/errors` -- which is what BayesTS 0.3.0 writes. A
  file written before that is migrated by running its forecast again. An object
  in a session is migrated the same way, by `add_posterior_forecasts()`, and
  says so rather than being read as one that was never forecast.

* **A folder of models is worked on without holding it.** `open_models()`
  returns a handle to a folder written with `write_to_hdf5()`, carrying a row
  per model -- where it is and what it is -- and none of the draws.
  `map_models()` reads one model, applies a function, writes it back and drops
  it, and `add_priors()`, `add_initial_values()`, `add_seed()`,
  `add_posterior_coefficients()`, `add_posterior_loglik()` and `thin()` have
  methods for the handle, so a step costs one model per worker rather than the
  whole folder. They take `models`, the models a step is applied to, which makes
  a long estimation resumable, and the seeds are the ones the same models get as
  a list, so a run taken in parts draws what a run over all of them draws. With
  `write = FALSE` `map_models()` only reads, which is how
  `selection_criteria()` compares a grid too large to hold, and
  `read_models_from_folder()` takes `draws` for reading part of every chain.
  This is the file-first workflow bgvars has for a global model, for any list of
  models: a lag and rank grid, or an expanding window.

* **BayesTS can be run on stored models.** `bayests_files()` returns a function
  that runs the BayesTS executable on a model file or on a whole directory of
  them and leaves the results where it wrote them. `bayests_posterior()` remains
  the way to draw a model that is in the session, but for a model that is
  already stored it copies the draws three times and holds them in R for no
  purpose; run on the files, nothing of them passes through the session. The
  draws are the same either way, since a model is drawn with the seed in its
  file.

* **A model can be analysed from its file.** `open_model()` returns a handle to
  a model written with `write_to_hdf5()`, carrying its specification, its data
  and the length of its chain but none of its draws. `map_draws()` reads the
  chain in pieces and applies a function to each, and `irf()` and `fevd()` have
  methods for such a handle: an impulse response stacks the responses of the
  pieces and takes its quantiles over all of them, a variance decomposition
  averages the pieces before shares are normalised or groups collapsed. Both
  give what the same call on the model in memory gives, which is what the tests
  check. This is for the models whose draws are what they are large in -- a time
  varying model, or a global model solved from many sub-models.

* **`read_model_from_hdf5()` reads part of a chain.** It read every draw of
  every block, so a caller working through a long chain -- solving a global
  model of time varying sub-models draw by draw, say -- had to hold the whole
  posterior of every model to get at one draw of it. The new `draws` argument
  takes the positions of the draws to read, and only those rows are read from
  the file; `integer(0)` gives the model, its data and its priors without the
  draws. A partial read cannot describe the chain it came from, so its blocks
  are labelled as a chain of their own, while a full read keeps the labels the
  file carries.

* **`add_predictive_loglik()` takes VAR models.** It refused anything but a VEC
  model, so the log predictive likelihood, the criterion for models whose
  coefficients or variances follow a state equation, was unavailable for the VAR
  models such a comparison puts a VEC model against. A VAR model is the case of a
  VEC model without an error correction term, so the density is the same
  expression with the cointegration block left out: windows of class 'bvarmodel'
  are now accepted, their rank counts as zero, and the check for a scaled or
  centred error correction term applies to VEC windows alone. Windows of mixed
  forms are still refused.

* **The vignettes and tests seed before `add_initial_values()`.** Since the
  seed of the posterior simulation is drawn by `add_initial_values()`, a
  `set.seed()` between it and `add_posterior_coefficients()` no longer changes
  the draws, and a model whose initial values were added first was simulated
  with a seed that no earlier call had fixed. The vignettes 'Introduction to
  bvartools', 'Bayesian Error Correction Models with Priors on the
  Cointegration Space', 'Bayesian Quantile VARs in bvartools' and 'Sign
  Restrictions in bvartools' and 26 places in the tests called `set.seed()`
  there.
  It now comes first, and the four vignettes are pre-compiled again.

* **Lists of models from other packages are seeded and simulated on several
  cores.** `add_seed()` on a `modellist` or an `expandingwindow` numbered only
  VAR and VEC models, and returned the dynamic factor models of dfmtools, which
  registers `add_seed()` methods of its own, without a seed. It now seeds and
  counts every element with an `add_seed()` method of its own, external
  forecasts excepted, and walks the same lists the parallel simulation does.
  With `cores` above 1, the workers loaded bvartools but never the package that
  registered the methods for such a model, so a list of dfmtools models stopped
  with "no applicable method". Each worker now loads, after the library paths
  of the session, the namespace of every method the models of the list
  dispatch to, and a model without a seed is given one through `add_seed()`
  whenever it has a method for it. `utils` is a new import.

* **Lists of models can be simulated on several cores.** The 'modellist' and
  'expandingwindow' methods of `add_posterior_coefficients()`,
  `add_posterior_forecasts()` and `add_posterior_loglik()` have a new argument
  `cores`, 1 by default, which keeps the sequential simulation. With more, the
  models -- counted through nested lists, such as the expanding windows of
  several specifications -- are shared out over a PSOCK cluster of that many
  workers, which is stopped when the call returns. The workers start with
  `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and
  `VECLIB_MAXIMUM_THREADS` set to 1, so that an optimised BLAS does not start
  one thread per core in each of them, and the session gets its values back
  afterwards. They use the library paths of the session. Coefficient draws are
  seeded per model and do not depend on the number of workers; a model without
  a seed is given one first. Forecasts draw from worker streams set up with
  `clusterSetRNGStream()` from R's generator, so `set.seed()` reproduces them
  for a given number of workers. An error on a worker is raised as on one core.
  For 16 TVP models over eight expanding windows, eight workers cut the time of
  `add_posterior_coefficients()` from 21 to 5 seconds. With OpenBLAS the
  parallel draws equal those of a session running its BLAS on one thread:
  OpenBLAS rounds some results differently on one thread than on several, and
  a Markov chain carries the difference forward. `parallel` is a new import.

* **The README explains how R's BLAS library affects speed.** A new section
  under "Installation" shows how to switch R to OpenBLAS on Windows and Linux,
  or to Accelerate on macOS. For a time varying parameter VAR this halved the
  sampling time and matched the standalone BayesTS executable, while compiling
  with `-O3` made no difference.

* **Every model carries the seed of its posterior simulation.**
  `add_initial_values()` stores a seed, drawn from R's random number generator,
  as `object$model$seed` unless the model has one, and the new `add_seed()`
  replaces it -- for a `modellist` or an `expandingwindow` with `seed`,
  `seed + 1`, ..., one per model, counted through nested lists. The internal
  samplers of `add_posterior_coefficients()` draw with that seed: R's generator
  is set to it, with R's default kinds, for the simulation and put back as it
  was afterwards. **This changes the draws of existing code**, not their
  distribution: a `set.seed()` between `add_initial_values()` and
  `add_posterior_coefficients()` no longer has an effect, and results of
  earlier versions do not reproduce draw for draw. `set.seed()` before
  `add_initial_values()` still makes a script reproducible, and a model now
  draws the same whichever worker of a cluster simulates it. A model without a
  seed draws from R's generator as before. `write_to_hdf5()` writes the seed as
  attribute `seed` of `/model` and `read_model_from_hdf5()` reads it back.
  `use_expanding_window()` counts a seed the model already has up from window
  to window.

* **Posterior simulation with the BayesTS executable.** `bayests_posterior()`
  returns a function for argument `posterior_function` of
  `add_posterior_coefficients()`. It writes the model to a file, lets the
  standalone `bayests posterior` draw into it with the model's seed and reads
  the draws back; nothing else about the model changes. The draws have the
  elements, dimensions and thinning of those of the internal samplers, whose
  C++ code BayesTS shares. A BayesTS build linked against OpenBLAS drew 3.4 to
  4.8 times faster than the internal samplers, which use R's reference BLAS,
  on constant, stochastic volatility and time varying VAR and VEC models of
  bgvars' `gvar2023` data.

* **`add_priors()` rejects non-positive Wishart degrees of freedom.** The
  documentation of `sigma$df` said "a non-negative integer" and `add_priors()`
  only stopped for negative values, but every Wishart sampler rejects
  `df <= 0`, so `sigma = list(df = 0, ...)` passed `add_priors()` and failed
  later in `add_posterior_coefficients()`. For VEC models this was hidden until
  the error prior was stored as given, since `add_priors()` used to add the
  rank to `df`. `add_priors()` now stops for `df <= 0` -- for a VAR after
  truncating `df` to an integer, as it is stored -- and the documentation says
  "a positive integer". The gamma `shape` keeps its non-negative check.

* **The error correction term can be centred.** `scale_error_correction()` has
  new arguments `scale`, `TRUE` by default as before, and `centre`, `FALSE` by
  default. `centre = TRUE` subtracts the sample mean of each stochastic series
  in the error correction term, leaves the restricted deterministic terms as
  they are and stores the means as attribute `"centre"`; with `scale = FALSE`
  the term is centred only. The unrestricted constant takes up the means, so
  centring needs one, and the starting values of the constant are shifted with
  the series. `rescale_error_correction()` shifts the draws of the constant back
  by `-alpha beta' m` -- period by period, with the loadings and cointegration
  vectors of that period, for a time varying model -- so the fitted values and
  the log-likelihood of every draw are unchanged, and gives the series their
  means back. For a model with time varying cointegration vectors this removes
  the part of a step of those vectors that acts as a random walk intercept on
  series far from zero. It does not settle how far the cointegration space
  moves: in the US sub-model of bgvars' `gvar2023`, centring alone roughly
  halved how far the error correction term and the constant moved over the
  sample and left the oil price equation at about 0.8 of the maximum likelihood
  residual standard deviation. `vec_to_var()` and `add_predictive_loglik()` refuse a centred term
  as they refuse a scaled one, `bvec()` adds the means back when it rebuilds the
  levels, and the means travel with a model written to HDF5.

* **The documentation says what decides the drift of a time varying
  cointegration space.** `?cointspace_prior`, the `coef$rate_alpha` item of
  `?add_priors.bvecmodel` and the TVP-SV VEC vignette put the residuals a time
  varying VEC model can absorb down to the drift of its loadings, and a small
  `coef$rate` or `coef$rate_alpha` as the remedy. The cointegration vectors take
  steps of unit variance that no rate reaches, and how far those steps move
  `Pi_t` is set by the prior precision of the loadings relative to the scale of
  `Pi` and by the levels of the series in the error correction term: for log
  levels far from zero a step acts as a random walk intercept, or the loadings
  are drawn close to zero and switch the term off. With the state variances of
  all coefficients pinned at 1e-14, the US sub-model of bgvars' `gvar2023`
  data set still reported the residual standard deviation of one equation as low
  as 0.1 to 0.5 of maximum likelihood in some chains, with which equation and
  how low depending on the chain. The section 'Prior on the cointegration space' now describes both
  channels, with the evidence, why scaling the error correction term makes them
  worse, and that the comparison of the residual variances with least squares
  remains the check to make.

* **Time varying cointegration spaces start on the scale of their state
  equation.** `add_initial_values()` started the cointegration vectors of a time
  varying VEC model at the maximum likelihood estimate normalised to
  beta' S11 beta = I, whose scale depends on the data and has nothing to do
  with the state equation beta_t = rho beta_{t-1} + eta_t, eta_t ~ N(0, I). Each
  vector now starts at the same direction rescaled to the stationary norm
  sqrt(k_w / (1 - rho^2)), and the loadings are estimated against it, so the
  starting Pi is unchanged; the starting values drawn from the prior get the same
  scale. On the Austrian sub-model of bgvars' TVP-SV-GVEC vignette the old start
  had a norm of 552 against a posterior of about 45, and the chain took some
  1,200 to 1,500 iterations to get there -- more than the vignette's burn-in of
  1,000 -- during which the cointegration space turned about ten times more
  slowly than it does in the posterior. From the new start, 59, it is there from
  the first draws. **Results change** for time varying VEC models, mostly in
  their early draws.

* **The loadings of a time varying VEC model can drift more slowly than its
  other coefficients.** New element `coef$rate_alpha` of `add_priors()` sets the
  rate of the gamma prior on the precisions of the loadings' random walks, as
  `coef$rate_det` does for the deterministic terms; without it the loadings take
  `coef$rate`, as before. The loadings multiply the levels in the error
  correction term, so one rate for all coefficients either freezes the
  coefficients of the differenced regressors or lets the loadings -- and, with a
  large `rate_det`, the constant -- absorb the residuals. On the sub-models of
  the TVP-SV-GVEC vignette of bgvars, a rate of 1e-5 with `rate_alpha = 1e-10`
  and `rate_det = 1e-8` let the short-run coefficients drift while every
  equation kept a residual standard deviation of the order of least squares,
  where a common rate of 1e-5 took it down to 0.3 to 0.8 of it. A VAR model
  refuses the element.

* **VEC models are forecast directly, with their states simulated forward.**
  `add_forecast_input()`, `add_posterior_forecasts()`, `predict()` and
  `add_forecast_errors()` now have methods for a `bvecmodel`; the first three
  used to stop and point to `vec_to_var()`. The forecasts are of the levels, as
  before: the regressors of the forecast periods are those of the VAR in levels,
  assembled from the data `vec_to_var()` builds, and `predict()` shows them with
  the history in levels. The simulation, though, runs on the VEC model's own
  draws through the vendored BayesTS core. For a model with constant coefficients
  that gives the same forecasts as converting first, from the same seed. For a
  time varying model it does not: `add_posterior_forecasts(forecast_states =
  "simulate")`, the default, steps the loadings, short-run coefficients,
  cointegration vectors, covariance block and log-volatilities of each draw and
  rebuilds the VAR in levels at every forecast period, where the VAR
  representation can only hold its coefficients at the last period.
  `forecast_states = "hold"` keeps them there. An expanding window of VEC models
  can now be forecast and evaluated with `selection_criteria()` without
  converting its windows. `irf()`, `fevd()` and `spillover()` still take the VAR
  representation.

* **Forecasts of VEC models with stochastic volatility carry the volatility
  forward.** `vec_to_var()` marked the VAR representation of every VEC model with
  stochastic volatility as `forecast_states = "hold"`, so its forecasts kept the
  volatility of the last period for every forecast period. With constant
  coefficients that representation is a `VarNormalStochvol` model with the same
  covariance block and the same log-volatility random walk, and its posterior
  carries the innovation variances, so it now simulates the volatility forward
  like any other VAR model with stochastic volatility. A model estimated before
  `posterior$u_sigma_inv$sigma` was kept still holds, and so does every time
  varying VEC model, whose level coefficients are not a random walk of their own.
  **Results change** for forecasts of VEC models with stochastic volatility and
  constant coefficients obtained through `vec_to_var()`; `add_posterior_forecasts(
  forecast_states = "hold")` gives the old ones.

* **Vendored BayesTS core refreshed: VEC samplers can simulate their states
  forward.** The core's forecasts of time varying and stochastic volatility VEC
  models now step the loadings, short-run coefficients, cointegration vectors,
  covariance block and log-volatilities and rebuild the level VAR at every
  forecast period, as its VAR forecasts already did. Nothing changes at the R
  level: a `bvecmodel` is still forecast through `vec_to_var()`, whose objects
  hold their states. The internal VEC bindings hand the core the draws a
  simulated forecast reads, so they stay usable from inside the package.
  Coefficient draws, log likelihoods and every forecast this package produces
  are unchanged.

* **Forecasts of time varying and stochastic volatility models carry the drift
  forward.** `add_posterior_forecasts()` forecast `VarTvpGamma`, `VarTvpStochvol`,
  `VarTvpWishart` and `VarNormalStochvol` models from each draw's coefficients,
  covariance block and volatilities in the last sample period, held for every
  forecast period. That is the forecast of a model whose drift stops where the
  sample does: its intervals left out the drift, and a held volatility
  understated the expected variance of every period after the first. The
  vendored BayesTS core now takes one step of each draw's random walks per
  forecast period, by the innovation variances the sampler drew for them, and a
  coefficient that BVS excluded stays at zero. The new argument
  `add_posterior_forecasts(forecast_states = "hold")` gives the old forecasts and
  is stored in `model$forecast_states`. For stochastic volatility models
  `add_posterior_coefficients()` now also keeps `posterior$u_sigma_inv$sigma`,
  the variance of the log-volatility innovations, which the forecast steps the
  volatility by; a model fitted with an earlier version lacks it and forecasts
  only with `forecast_states = "hold"` until it is fitted again. **Results
  change** for the forecasts of those four models. Their coefficient draws and
  log likelihoods are unchanged, which upstream verified by fingerprinting every
  sampler before and after. VEC models still hold their states when forecast,
  but `VecNormalStochvol` and `VecTvpStochvol` keep `posterior$u_sigma_inv$sigma`
  as well, so their posteriors carry the same blocks as those of the VAR models.

* **Predictive likelihoods for VEC models.** New `add_predictive_loglik()` adds
  to every window of an expanding window exercise the draws of the one-step-ahead
  log predictive density of the observation the next window adds, and
  `selection_criteria()` sums them to the log predictive likelihood, criterion
  `"LPL"`, which `choose_best_model()` maximises. Time varying coefficients,
  cointegration spaces, error covariances and stochastic volatilities are
  carried one period forward with their state equations first. The variance of
  the random walk of the log volatilities is read from
  `posterior$u_sigma_inv$sigma` where a sampler stores it, and drawn from its
  conditional posterior given the path otherwise. This is the criterion Koop,
  León-González and Strachan (2011) choose between time varying cointegration
  models and ranks with: the pointwise log-likelihood behind WAIC and LOOIC is
  evaluated at states that have seen the period's observation, and on the
  Austrian sub-model of a global VEC model it left the ranks of a TVP-SV-VEC
  0.2 apart where the predictive likelihood separated them by 3.7 with a standard
  error of 1.5. `?selection_criteria` says when to use which.

* **The prior of the loadings of a time varying VEC model follows `coef$v_i`.**
  `add_priors()` gave the loadings a prior precision of 1 / (1 - rho^2)
  whatever `coef$v_i` was, which put alpha beta' at unit prior variance rather
  than on the scale `coef$v_i` asks for, as `?cointspace_prior` says. Koop et al.
  (2011) scale the variance of the loadings by 1 - rho^2 relative to that of
  the other coefficients, and the precision is now `coef$v_i / (1 - rho^2)`.
  **Results change** for time varying VEC models with a `coef$v_i` other than 1;
  the tvp-sv-vec vignette uses 1 and is unaffected.

* **Time varying coefficients leave their starting values.** The time varying
  samplers of the vendored BayesTS core drew a coefficient path against the
  previous draw of the state before the sample, with the random walk's own
  innovation variance as the prior covariance of the first period, and then drew
  that state given the path. The two steps tied each other with that variance,
  so with a small one -- a `coef$rate` of 1e-8 or 1e-12, as the time varying
  vignettes use -- the chain could not move the level of a coefficient path:
  the posterior of the coefficients was the output of `add_initial_values()`.
  Two chains started from different values returned their own starting values,
  and a time varying VEC model kept its loadings at their maximum likelihood
  starting values while its cointegration vectors moved, which could leave it
  without error correction. The path is now drawn with the state before the
  sample integrated out of the first period's prior, and that state is drawn
  before the innovation variance. This covers the coefficients and the
  covariance block of `VarTvpGamma`, `VarTvpStochvol`, `VarTvpWishart`,
  `VarTvpAld`, `VecTvpGamma`, `VecTvpStochvol` and `VecTvpWishart`. The prior
  precision of the state before the sample must now be positive definite, since
  the draw takes its inverse, so `coef$v_i = 0` stops a time varying model with
  a message saying so. **Results change** for every time varying model.

* **Forecasts of an expanding window start from the end of each window.**
  `prepare_forecast_input()` took the lags of the first forecast period from
  the original series a model was created from rather than from its estimation
  sample, so every window of `use_expanding_window()` but the last forecast from
  the last observations of the whole series, and so did a model cut short by
  `window()`. On the first differences of `at_macrodata` the lagged values were
  off by up to 7.9 for a structural VARX, and by about 1 for a VEC model in
  levels through `vec_to_var()`. The forecast errors and out-of-sample criteria
  of expanding windows were computed from these forecasts. **Results change**
  for the forecasts of every window but the last.

* **`add_initial_values(method = "prior")` draws from the prior.** The
  coefficients, the covariance coefficients and the initial log-volatilities
  were drawn with the prior precision as their covariance, so a prior standard
  deviation of 0.1 started the chain at coefficients with a standard deviation
  of 10. The precisions of a gamma error prior were drawn as the inverse of a
  gamma with half the shape and the inverse rate, near 0.001 instead of near 2
  for `shape = 50` and `rate = 25`, and the state precisions of time varying
  models from a gamma with half the shape and rate, while the samplers read
  `shape` and `rate` as a Gamma(shape, rate) prior on a precision. A prior that
  is uninformative for some coefficients now stops with a message saying so,
  where `chol()` reported a leading minor that was not positive. Only the
  starting point of the chain changes.

* **Semiautomatic SSVS priors for structural models.** `add_priors()` refused
  `varsel$semiautomatic` for structural models, and `ssvs_prior()` fell back to
  `tau` for the contemporaneous coefficients of a structural VEC model. Both
  now scale by the least squares standard errors of the recursive system,
  equation by equation: the regressors of equation i together with minus the
  current values of the variables before it, and the error variance of each
  equation with its own degrees of freedom.

* **The Minnesota-like inclusion prior covers the contemporaneous coefficients
  of structural models.** They kept `inprior` whatever `varsel$minnesota` said,
  and now get kappa2, the inclusion probability of the other endogenous
  variables, as the Minnesota prior gives them the variance of the other
  endogenous variables. **Results change** for structural models with variable
  selection and `varsel$minnesota`.

* **Gamma priors are documented, and can differ by equation.** `?add_priors`
  now says that `shape` and `rate` are those of a gamma prior on a precision,
  with mean `shape / rate`, which is how every sampler reads them; a comment in
  the code claimed the error prior was halved on its way to the sampler, which
  it is not. `sigma$shape` and `sigma$rate` of a gamma or asymmetric Laplace
  error prior can be given per equation, which stopped with "the condition has
  length > 1".

* **`predict()` returns every simulated period by default.** `n_ahead`
  defaulted to 10 and warned whenever fewer periods had been simulated. It now
  defaults to the horizon given to `add_forecast_input()`.

* **`write_to_hdf5()` writes the starting error precision of constant gamma
  models where BayesTS reads it.** R keeps the starting precision of every
  model with a gamma error as `initial$u_omega_inv`, and the model file carried
  that name. BayesTS reads it from `/initial/u_sigma_inv` for `VarNormalGamma`
  and `VecNormalGamma`, and from `/initial/u_omega_inv` only for the time
  varying gamma models, so it refused every exported VAR or VEC with constant
  coefficients and an `error` of `"gamma"` or `"gamma+covar"`, structural models
  among them, with "initial error precision must be 6x6, got 0x0". The file now
  uses the BayesTS name for these models, and `read_model_from_hdf5()` maps it
  back, so a model read from a file still estimates inside R. Models estimated
  inside R were not affected.

* **`summary()` prints models with a single endogenous variable.** Printing
  selected the variance from the lower triangle of the covariance table, which
  for one equation is a single row and was dropped to a vector, so the credible
  bounds were looked up as columns that no longer existed and printing stopped
  with "undefined columns selected". Every AR model and every single-equation
  quantile regression failed this way. The same selection in the summary of VEC
  models is fixed too.

* **`add_priors()` stores the error prior of VEC models as given.** It added
  the cointegration rank `r` to `sigma$df` and to the gamma `sigma$shape`. The
  VEC samplers with constant coefficients and a Wishart prior add `r` to the
  posterior degrees of freedom themselves, which the prior of the loadings
  given the error covariance requires (Koop, León-González and Strachan, 2010,
  eq. 8), so `r` entered twice, and for the gamma shape and the time varying
  models there was no reason for it. `sigma$df` and `sigma$shape` now mean the
  same for VEC models as for VAR models. **Results change** slightly for every
  VEC model with a positive rank and a Wishart or gamma error prior: its prior
  on the errors is `r` degrees of freedom weaker.

* **VEC models with restricted deterministic terms or exogenous variables sample
  the posterior of their cointegration space prior.** The samplers of the VEC
  models with constant coefficients follow Koop, León-González and Strachan
  (2010), who derive them for a cointegration term with one row per endogenous
  variable. A constant, trend or seasonal dummies restricted to the cointegration
  space, or exogenous variables, add rows, and the draw of the cointegration
  matrix then left out a factor of the prior. The posterior overstated the size
  of `Pi`: in a simulation-based calibration with two endogenous variables and a
  restricted constant, the true size of `Pi` had a mean rank of 0.39 among its
  posterior draws over 1000 replications, where a correct sampler gives 0.50. The
  vendored BayesTS core now draws it exactly, by giving the loadings the rows
  they lack for the duration of that draw, after which the calibration gives mean
  ranks of 0.49 to 0.51 with restricted constants, trends and exogenous
  variables. The draw is a Gibbs step, so the cointegration matrix moves in every
  iteration: in a VEC model of Austrian output, inflation and short and
  long-term rates from `at_macrodata`, with three foreign series, a restricted
  trend and a flat prior, the posterior mean of `Pi` lies within 0.8 posterior
  standard deviations of its maximum likelihood estimate, and the median
  effective sample size of its elements is about 900 of 4000 draws.
  **Results change** for VEC models created with
  `const`, `trend` or `seasonal` set to `"restricted"` or with `exogen`, unless
  `tvp = TRUE`. Models without them are unaffected, and their draws are the same
  as before for a given seed.

* **Minnesota-like inclusion priors work for VEC models with exogenous
  variables.** `inclusion_prior()` with `minnesota_like = TRUE`, and so
  `add_priors()` with `varsel$minnesota`, filled `s` blocks of lagged exogenous
  differences after the current one, although a VEC model has `s - 1`. Whenever
  the model had fewer unrestricted deterministic terms than exogenous variables,
  as with the three foreign series of `at_macrodata` and a constant, this
  stopped with "subscript out of bounds". With at least as many deterministic
  terms the surplus block was overwritten by theirs, so the probabilities were
  already right there.

* **`add_initial_values()` stores the LS error covariance coefficients in the
  order the samplers read them.** The strict lower triangle of `Psi` is stored
  row by row, but the regression that estimates it returned its coefficients
  column by column. The two orders agree for up to three endogenous variables
  and differ from four on, so a model with `error = "gamma+covar"` or
  `"sv+covar"` started its chain with the coefficients in the wrong places.
  A constant coefficient model soon draws its way out of that. A time varying
  one whose state variances start small, as `add_initial_values()` sets them
  under a small prior rate, stays at those starting values, and so did its
  error variances: in a six-variable `sv+covar` model the error standard
  deviations of the last three equations came out 9 to 52 times their least
  squares values, and its log-likelihood, WAIC and LOOIC with them. **Results
  change** for models with a covariance block, at least four endogenous
  variables and the default `method`.

* **VEC models with constant coefficients no longer stop with
  `sqrtmat_sympd(): transformation failed`.** The vendored BayesTS core
  normalises the cointegration draws through a singular value decomposition
  instead of the matrix square root of their cross product, which failed on an
  ill-conditioned draw: a full-rank VEC model with a Wishart prior and BVS on
  six Austrian series stopped partway through its chain. Draws change by a
  rounding error only, at most 1.2e-13 relative in BayesTS's test suite.

* **`create_bvecmodel()` requires the lag orders of the VAR in levels.**
  Argument `p` no longer defaults to 2, and `s` must be given whenever
  `exogen` is. Both are lag orders of the corresponding VAR in levels, as is
  common in the literature, so the VEC model has `p - 1` lags of differences
  and the exogenous variables enter with their contemporaneous difference and
  `s - 1` lags of it. Leaving them without a default avoids a VEC model whose
  lags were never chosen, and avoids reading `p` as the number of lagged
  differences. `p` and `s` must be whole numbers of at least 1, and `s` is
  ignored without `exogen`. **Breaking change:** calls that relied on the
  defaults now stop with an error.

* **`minnesota_prior()` gives exogenous variables of VEC models their own prior
  variances.** The method for `bvecmodel` counted one block of exogenous
  differences less than the model has before it placed the deterministic
  terms, so their prior variance, `kappa1 * kappa4` times the residual
  variance, overwrote the last block of the exogenous variables. With `s = 1`
  that was every exogenous coefficient: with `us_macrodata`, `Dp` and `u`
  endogenous and `r` exogenous, the prior variance of the difference of `r`
  in the equation of `Dp` came to 1.76, the value of the deterministic terms,
  whatever `kappa3` said, and with `s = 2` the same happened to its first
  lag. The regressions that give the residual standard
  deviations used only the first of several deterministic terms. **Results
  change** for Minnesota priors of VEC models with exogenous variables, and of
  VEC models with more than one restricted or unrestricted deterministic term.

* **`generate_artificial_var()` simulates time varying, stochastic volatility
  and structural models.** With `tvp = TRUE` the coefficients and the free
  elements of Psi follow random walks, with `sv = TRUE` the log-volatilities
  do, and with `structural = TRUE` the series come from a lower triangular A0
  with uncorrelated structural errors. These are the state equations and the
  structural form that `create_bvarmodel()` estimates. Time varying parameters
  are returned as `K x M x T` arrays, whose `as.vector()` is in the order of the
  posterior draws, together with the variances of their state equations.
  Coefficients that are zero stay zero, and with `stable = TRUE`, the default,
  the process is stable in every period. The function also gains a linear
  trend, a lag order of zero and `stable = FALSE` for integrated or explosive
  series. **Breaking changes:** argument `const` is replaced by
  `deterministic`, which takes the values of `create_bvarmodel()`, and
  `a_range` is renamed to `range_a`. `range_const` now has its documented
  default `c(-0.5, 0.5)` instead of stopping when it is not given. The default
  `presample` is 100 instead of 0, so the series no longer depend on their
  starting values, and coefficients are no longer rounded to two digits.
  **Results change** for a given seed.

* **New function `generate_artificial_vec()`** simulates cointegrated series
  from a VEC model with rank `r`, restricted or unrestricted constants and
  trends, and the same options `structural`, `tvp` and `sv`. Its coefficients
  are drawn so that the VAR representation has exactly `k - r` unit roots, in
  every period for time varying models. Besides `alpha`, `beta` and `gamma` it
  returns `pi`, which estimates can be compared with whatever normalisation of
  `beta` they use.

* **Structural impulse responses and variance decompositions read A0 correctly
  for models with four or more endogenous variables.** The free elements of
  A0 are stored column by column, which is also the order of the
  contemporaneous regressors, but the function that prepares the draws for
  `irf()` and `fevd()` read them row by row. For up to three variables the two
  orders coincide. From four on, elements were swapped: with `us_macrodata` and
  the change in its interest rate, the six-quarter decomposition of that change
  came to 24, 21, 55 and 0 percent instead of 12, 22, 66 and 0 percent, and a
  structural generalised response had the wrong sign. Affected are types
  `"sir"` and `"sgir"`, for time varying models as well and for structural VEC
  models through `vec_to_var()`. **Results change** for structural models with
  four or more endogenous variables.

* **`summary()` and `plot()` show each element of A0 in its own cell.** They
  had the same row by row reading, so with four or more variables A0 and the
  inclusion probabilities of its elements were reported in each other's cells.

* **Structural VEC models can be summarised.** `summary()` stopped with
  "update structural" for every structural VEC model, a placeholder left in the
  code. Its structural branch now counts the regressors of the VEC model, names
  the columns of A0 after the endogenous variables, and removes the loadings
  from the inclusion probabilities before it places those of A0.

* **`fevd(type = "gir")` divides each shock by its own variance, as in Pesaran
  and Shin (1998).** The generalised decomposition divided every share by the
  standard deviation of the *response* instead, which is the decomposition of
  no shock at all. Normalising did not repair it, because the scaling that was
  missing differs from shock to shock. In a VAR of output growth, inflation and
  the short-term interest rate from `at_macrodata`, the eight-quarter shares of
  inflation came to 46, 53 and 1 percent with `normalise_gir = TRUE`, where
  Pesaran and Shin give 9, 77 and 14 percent. `type = "sgir"` had the same
  scaling and is corrected with it. The documentation, which gave a third
  version dividing by the variance of the response, now states the formula the
  code uses. The shares agree with `spillover()`, which already used the
  correct scaling. **Results change** for every generalised decomposition.

* **Orthogonalised and generalised impulse responses are responses to shocks of
  one standard deviation.** With the default `shock = 1`, `irf(type = "oir")`
  scaled the Choleski factor to a unit diagonal and `type = "gir"` divided by
  the variance of the impulse variable, so both returned responses to a unit
  shock. The documentation gave the one standard deviation versions,
  `Phi_i P` and `sigma_jj^(-1/2) Phi_i Sigma e_j`, which is also what
  `vars::irf()` returns. The code now follows the documentation, and
  `shock` counts standard deviations for `"oir"`, `"gir"` and `"sgir"`, so
  `shock = "sd"` and `"nsd"` are the same as `1` and `-1` there. `"feir"` and
  `"sir"` still count units of the error. Orthogonalised responses are now on
  the scale of sign restricted ones, whose impact matrix is a rotation of the
  unscaled Choleski factor. **Results change** for `"oir"`, `"gir"` and
  `"sgir"` responses with a numeric `shock`.

* **AIC, BIC and HQ start from the deviance at the posterior mean.** The
  deviance at the point estimate was recovered from the draws as their mean
  deviance less the sum of the pointwise variances of the log-likelihood, the
  penalty of WAIC. That sum is not the effective number of parameters. For VARs
  of lag orders 0 to 4 with a flat prior on `at_macrodata`, which have 9 to 45
  parameters, it came to 11 to 58, so AIC fell 2 to 11 below its maximum
  likelihood value, by more the larger the model, and AIC ranked the lag orders
  3 and 4 the wrong way round. The log-likelihood of the model is now evaluated
  once at the posterior mean of its parameters. For a VEC model the point is
  the matrix of the rank of the model closest to the posterior mean of
  `Pi = alpha beta'`, because `alpha` and `beta` are identified only up to a
  rotation. Closest is measured as the likelihood measures it, weighting a
  difference in `Pi` with the error correction term and the error precision,
  so the point does not depend on the scale of the series. Measuring every
  element of `Pi` alike, as a first version of this change did, gave up fit
  where the levels are large: for VEC models of four `at_macrodata` levels of
  rank 0 to 3, AIC of rank 2 came to 2213 against -121 at the maximum
  likelihood estimate, and the ranks were ordered unlike maximum likelihood.
  With the weighted point AIC lies within 7 of its maximum likelihood value and
  orders them alike. The note about periods with a highly variable pointwise
  log-likelihood now refers to WAIC alone. **Results change** for AIC, BIC and
  HQ of every model.

* **`spillover()` decomposes the `n_ahead` step forecast error variance.** Its
  documentation and Diebold and Yilmaz (2012) sum the impulse responses of
  periods 0 to `n_ahead - 1`. The code summed them up to `n_ahead`, so the
  default of 10 decomposed the 11 step error. The table now corresponds to
  period `n_ahead - 1` of `fevd()`, which counts its periods from impact, and
  the documentation of `fevd()` states that its period `h` sums the
  responses of periods 0 to `h`. **Results change** for every spillover
  measure.

* **Forecast input for a model with exogenous variables says which periods
  `exogen` has to cover.** The regressors of a forecast period include the lags
  of the exogenous variables, so `exogen` has to reach back `s` periods into
  the estimation sample. Nothing documented this -- the argument pointed to a
  'Details' section that did not exist -- and a series starting with the first
  forecast period failed with R's "number of items to replace is not a multiple
  of replacement length". `add_forecast_input()` now refuses such a series with
  a message naming the periods it has to cover, and `?prepare_forecast_input`
  explains why.

* **Time varying models repeat their posterior draws from a seed set after
  `add_initial_values()`.** For every model with `tvp = TRUE`,
  `add_initial_values()` drew the initial precisions of the state equations,
  `a_sigma_inv` and `psi_sigma_inv`, from their gamma priors under every
  `method`. Under the default `"ols"` or `"maxlik"`, that was its only use of
  the random number generator. A `set.seed()` placed between it and
  `add_posterior_coefficients()` fixed the sampler but not where the chain
  started, so two R sessions running the same script gave different draws.
  This looked like undefined behaviour in the samplers, most visibly for
  `VarTvpAld` and `VarTvpGamma` with `varsel = "bvs"`, but the samplers
  themselves repeat exactly given the same inputs. The precisions now start
  at the prior mean, `shape / rate`, and are drawn only under
  `method = "prior"`. **Initial values and draws change** for every time
  varying model estimated with `method = "ols"` or `"maxlik"`.

* **Vendored BayesTS core refreshed again: the covariance block of the time
  varying gamma models, BVS in the quantile VAR, and validation of prior
  values.** Upstream fixed what four further audits found. The same thirty-six
  specifications as for the previous refresh, now including BVS on every time
  varying and quantile model, were fitted from pinned seeds against the package
  built before and after, and compared block by block.

    **`VarTvpGamma` and `VecTvpGamma` draw the covariance block under the
    current error variances.** They inverted the starting error precision once,
    before the chain, and drew every path of the covariance block under that
    inverse, although the precision was redrawn in every iteration. **Draws
    change** for both models with `error = "gamma+covar"`, with and without
    variable selection.

    **BVS in `VarNormalAld` uses the data to bring an excluded coefficient back.**
    Its likelihood masked each candidate a second time, with the indicators the
    sweep was still updating, so a coefficient that was out came back on the
    prior inclusion probability alone. **Draws change** for
    `create_bvarmodel(error = "ald", varsel = "bvs")`.

    **`VarNormalStochvol` no longer fails on an observation far in the tails.**
    It carried its own copy of the mixture draw of the log volatilities, whose
    component probabilities could all underflow to zero and whose index could run
    off the end of the table. It now uses the shared draw that the other
    stochastic volatility samplers use. The model is the same; here the draws
    differ from before only by rounding (at most 2e-9 relative).

    **The samplers refuse prior values no model can mean**: a negative or
    non-finite gamma shape or rate, an inclusion probability outside [0, 1], an
    asymmetric prior precision or Wishart scale, and a non-diagonal starting
    precision where only the diagonal is redrawn. VAR forecasts also refuse
    regressors without exactly one row per horizon, where five of them used to
    read past the end.

    **VEC models with BVS or SSVS can be estimated again.** The new check on
    inclusion probabilities refused every one of them, because
    `inclusion_prior()` set the prior inclusion probability of the loadings to
    `NA`. The loadings are never selected, so the value was never read. It is now
    1, which says what happens to them. A test now fits a VEC model with each of
    the two selection methods, which no test did before. The draws of these
    models are unchanged.

    Every other specification is bit-identical, forecasts included, with two
    exceptions this refresh did not cause. A time varying quantile VAR with BVS
    did not repeat its draws from the same seed in a fresh session, and a time
    varying VAR with a gamma error term and BVS differed between the two builds.
    Neither was the samplers: `add_initial_values()` drew the starting state
    precisions without a seed, which the entry above fixes.
    The vendored set is still the same 65 files.

* **Vendored BayesTS core refreshed: BVS, SSVS and the log likelihood of the
  time varying models are fixed.** Upstream found three errors in an audit, and
  all three reached the models of this package. Thirty specifications -- every
  VAR and VEC algorithm, each variable selection it offers, with and without a
  covariance block -- were fitted from pinned seeds against the package built
  before and after the refresh, and compared block by block. What moved is
  listed below. Nothing else did.

    **BVS draws its inclusion indicators from their posterior.** The sweep set an
    indicator to one with probability `min(1, exp(l1 - l0))` rather than the
    logistic of `l1 - l0`. It also scored every candidate against the inclusion
    matrix as it stood before the sweep. The rule came from this package's
    original `bvs.cpp`, and the chain it made did not have the posterior as its
    stationary distribution: at a prior inclusion probability of 0.5 and an
    uninformative likelihood, a coefficient once included stayed in. Each
    indicator is now a Gibbs draw given the current state of all the others
    (Korobilis, 2013). Separately, `VarNormalGamma`, `VarNormalStochvol`,
    `VecNormalGamma` and `VecNormalStochvol` applied BVS to the covariance block
    against regressors copied before the chain started, while they were still
    zero, so `varsel$covar = TRUE` selected on the prior alone. **Draws change**
    for every model with `varsel = "bvs"`, in every block.

    **SSVS includes a coefficient far from zero.** The spike and slab densities
    were formed before they were divided, so for a coefficient many slab widths
    from zero both were zero, the inclusion probability was `NaN`, and the
    coefficient the data most want in was excluded. The odds are now formed in
    logs. **Draws change** for every model with `varsel = "ssvs"`.

    **The time varying models score every period of the log likelihood under
    its own error precision.** `add_posterior_loglik()` scored the whole sample
    under the precision of the last period, which is the likelihood of a model
    whose volatility does not move. That was wrong for `error = "sv"` and
    `"sv+covar"`, and for `"gamma+covar"`, where the covariance block makes the
    precision drift. WAIC and LOO of such a model, and therefore the model
    comparisons built on them, compared the wrong thing. Taking the fix needed a
    change here as well: the log likelihood bindings of `VarTvpGamma`,
    `VarTvpStochvol`, `VecTvpGamma` and `VecTvpStochvol` in `src/` sliced the
    precision to the last period before handing it over, and now pass the whole
    path. Their forecasts still start from the last period. **Draws are
    unchanged; the log likelihood changes** for those models.

    The determinant term of every Gaussian log likelihood now goes through a
    Cholesky factor in logs, so it cannot underflow for a large model. Elsewhere
    the log likelihood changes by at most 3e-11 relative, and coefficients and
    forecasts are bit-identical. The constant VEC samplers also refuse a prior
    precision of the cointegration space that is not symmetric, which
    `add_priors()` already did.

    The refresh adds one file, `src/core/algorithms/inclusion_probability.h`,
    which `inst/COPYRIGHTS` now lists.

* **The package has a DOI, 10.5281/zenodo.22736604.** The GitHub release of
  0.3.0 is archived on Zenodo, and this concept DOI identifies the package as a
  whole: <https://doi.org/10.5281/zenodo.22736604> resolves to whichever version
  was archived last, so it stays valid from one release to the next.
  `citation("bvartools")` prints and exports it, `CITATION.cff` carries it for
  the "Cite this repository" box on GitHub, the README shows it as a badge,
  and every vignette ends with a short section on citing the package.

* **The documentation warns against scaling the error correction term of a VEC
  model with time varying coefficients.** In the vignette on TVP-SV-VEC models,
  full-length chains on a scaled error correction term let the coefficient paths
  absorb the residuals of the short-term interest rate equation in one chain out
  of four, depending on nothing but the seed: the volatility of that equation
  collapsed to a nearly flat fraction of the one of a least squares fit. A tighter
  `coint$rho` made that rarer without preventing it. On unscaled series with
  `coef$rate = 0.000001`, eight chains with different seeds agreed. The section on
  the prior on the cointegration space in `?cointspace_prior`, which
  `?add_priors.bvecmodel` shows, now says so. It recommends unscaled series, a
  small `coef$rate` and a comparison of the residual variances with those of a
  least squares fit across seeds. The vignette estimates its models that way.

* **The prior on the cointegration space is built by an exported function,
  `cointspace_prior()`.** `add_priors()` for VEC models checks `coint` and builds
  `priors$beta` through it, and a package with VEC models of its own layout can
  call it too -- bgvars does for the sub-models of a global VEC model, whose error
  correction term also holds the weakly exogenous and global variables. They then
  get the same checks and the same prior, the one centred on the maximum likelihood
  estimate included, instead of a copy that could drift from it. Its dimensions come
  from the error correction term `data$train$w` rather than from `model$m`, which a
  model of another layout counts differently. For a sub-model with lagged foreign
  variables that gave a numeric `p_tau_i` a matrix of the wrong size. For the models
  of `create_bvecmodel()` the two agree, and the priors are unchanged bit for bit.
  The documentation of `coint` moved to `?cointspace_prior` and is shown in
  `?add_priors.bvecmodel` as a section of its own.

* **`add_posterior_forecasts()` and `add_posterior_loglik()` no longer change
  the object they are given, and a second forecast replaces the first.** The
  C++ functions behind them take the model as an `Rcpp::List`, which wraps the R
  object the caller holds rather than a copy of it, and they wrote their result
  into that list. After `add_posterior_loglik(object)` without an assignment,
  `object` had gained a `loglik` all the same. The result was also appended
  rather than set by name. Forecasting an object that already carried a
  forecast -- under a new seed, or after its coefficients were edited in R --
  therefore left two `forecast` elements, and `object$posterior$forecast`
  silently returned the first of them. The log likelihood behaved the same way.

    All thirteen VAR and VEC bindings that forecast and all fifteen that compute
    a log likelihood now return through `with_posterior_element()` in
    `src/bayests_r_io.h`. It copies the model and its posterior list before it
    sets the element by name. The copy is shallow, so the draws already stored
    are shared with the caller's object rather than duplicated. How the draws
    are computed is unchanged. The `*Coefficients` bindings already built a new
    list and needed no change.

* **A VEC model with constant coefficients and stochastic volatility can be
  forecast.** `vec_to_var()` stopped on such a model with "posterior draws of
  u_sigma_inv must have 9 rows, got 1584": the coefficients are constant, so it
  handed them to the transformation in one piece, and with them the whole path
  of the error precision, which the transformation checks against the size of
  a single period. It now passes the precision of the last period, which the
  coefficients do not depend on, and carries the path over to the VAR
  representation unchanged. Its forecasts start from the precision of the last
  in-sample period, as those of the time varying models do. VEC models with a
  constant error covariance and those with time varying coefficients were not
  affected.

* **The samplers thin: `create_bvarmodel()` and `create_bvecmodel()` take
  `thin`.** A model with `thin = t` runs its chain for `burnin + iterations * t`
  draws and keeps the last of every `t` after the burn-in. So `iterations` is
  still the number of draws kept, and the posterior and every file written from
  it are still sized by it. `thin()` still thins draws after they are made; the
  argument is what lets a slowly mixing chain run long without holding every
  draw in memory. A chain thinned this way is exactly every `t`-th draw of the
  unthinned chain from the same seed, and a test asserts that for a VAR and a
  VEC.

    The `mcpar` of the draws says which ones were kept: iterations `t`, `2t`,
    ... after the burn-in, which is also what `bayests` writes into a model file.
    The log likelihood and the forecasts take it from the coefficients as
    before. `thin()` on a thinned model now counts from those labels instead of
    restarting at one, so that thinning a `t`-thinned chain by 2 labels its
    draws `2t`, `4t`, .... `model$thin` is carried only when it is above 1, the
    way `quantile` is carried only by the quantile models, so a model created
    without it is unchanged, and `vec_to_var()` passes it on.

    **Draws are unchanged** for every model that does not thin: the sampler
    reads a `thin` of 1, and the draws come back as the same `mcmc` objects as
    before.

* **Vendored BayesTS core refreshed: the samplers can thin.** **Draws are
  unchanged**, verified here rather than taken on trust. Sixteen specifications
  -- all fifteen VAR and VEC algorithms the package reaches, plus
  `VarNormalGamma` in both its structural and its covariance block form -- were
  fitted from pinned seeds against the package built before and after the
  refresh. The comparison covered the coefficients, the log likelihood and a
  four-step forecast wherever the model has one. Every posterior block of every
  one of them is bit-identical.

    Upstream's `VarSpec` gained `thin`, with `keeps()` and `kept_index()`, and
    all of its samplers now keep their draws through those two calls instead of
    testing `draw >= burnin`. At `thin = 1` the two are exactly that test, which
    is why nothing moved. Nothing in this package sets `thin` yet: `read_spec()`
    in `src/bayests_r_io.h` does not read it, so every model still keeps every
    draw after the burn-in, and `thin()` still thins after the fact. Passing a
    thinning interval to the sampler, so that a long chain is never held in
    memory whole, is a change for another time.

    The vendored set is still the same 64 files, so `inst/COPYRIGHTS` is
    unchanged.

* **`write_to_hdf5()` releases the file of a model before it returns.** The
  specification of a model, its class and the properties of its datasets are
  stored as HDF5 attributes, and they were written with hdf5r's `h5attr<-`,
  which never closes the attribute it creates. HDF5 keeps a file open for as
  long as anything in it is, so closing the file did not release it, and the
  handles were left to R's garbage collector. Until it ran, the file was locked
  against every other process on Windows: `bayests` could not open a model that
  had just been exported and failed with "Unable to lock file", and succeeded
  on the unchanged file after `gc()`. Each attribute is now closed as soon as it
  is written. Reading a model was not affected.

* **`create_bvarmodel()` uses `s = 2` by default, as its documentation already
  said and as `create_bvecmodel()` does.** The signature had `s = NULL`, so a
  call with `exogen` but without `s` stopped with "attempt to set an attribute
  on NULL" after a warning from `max()`. Such a call now includes the external
  regressors with two lags. An `s` that is not a vector of non-negative integers
  is refused with a message naming the argument.

* **Smaller fixes from an audit of the package against the TVP-SV VEC
  vignette.**

  - `summary()` accepted a `period` of 1.5 and silently used period 1, and
    stopped with the message of an `if()` on a vector of periods. `period` must
    now be a single integer within the sample, and the same check applies to the
    draws behind `irf()` and `fevd()`. `thin()` stopped with "wrong sign in 'by'
    argument" for an interval longer than the chain; the interval must now be a
    positive integer no larger than the number of draws.
  - `irf()`, `fevd()`, `spillover()`, `add_posterior_forecasts()` and `predict()`
    on a `bvecmodel` stopped with "no applicable method". Like
    `add_forecast_input()`, they now say to use `vec_to_var()` first.
  - The draws of `rho`, where a time varying VEC model puts a prior on it, are
    summarised by `summary()`, and `print()` says whether `rho` is fixed or
    drawn, and from which prior.
  - A round trip through `write_to_hdf5()` and `read_model_from_hdf5()` now
    returns every element that was written, with its values, shape and class;
    only the order of the elements within a group follows the file, which lists
    them by name. Scalar and character priors came back
    as 1 x 1 matrices, `priors$u_sigma$type` was not written, the class of a
    series changed in whether it listed "array", `model$rclass` was added, and
    list elements holding `NULL` were dropped. Values without dimensions are
    marked in the file so that they are read back as vectors; files written
    elsewhere carry no mark and are read as before. The `NULL` placeholders are
    no longer created in the first place.

* **Time varying VEC models and expanding windows of VEC models can be
  analysed through their VAR representation.** `vec_to_var()` used to refuse
  the posterior draws of a model with time varying parameters. It now applies
  the transformation period by period and returns a VAR model with time
  varying parameters, whose forecasts start from the last period and whose
  impulse responses and variance decompositions are those of the period asked
  for. The draws of the state variances of the VEC coefficients and of `rho`
  are dropped. A new `vec_to_var()` method for expanding windows opens the
  route to forecasts, forecast errors and out-of-sample selection criteria for
  VEC models, which `add_forecast_input()` pointed to but no method provided.

* **Unreliable WAIC, AIC, BIC and HQ are flagged.** All four rest on the
  variance of the pointwise log-likelihood across draws -- WAIC as its penalty,
  AIC, BIC and HQ through the deviance at the point estimate. Where that
  variance exceeds 0.4 in a period the correction cannot be relied on (Vehtari,
  Gelman and Gabry, 2017), and a time varying or stochastic volatility model
  sampled briefly could report criteria far below its fit, even negative, with
  nothing to say so. The number of such periods is now stored with WAIC as
  attribute `n_high_var` and printed under the criteria, as the Pareto k
  diagnostic of LOOIC already was.

* **`add_priors()` rejects unknown elements in all of its arguments.** It
  already did for `coef`, but ignored a misspelt element of `sigma`, `varsel`
  and `coint` without a message, so the setting it was meant for was silently
  left out. The horse-races vignette, for instance, passed
  `varsel$exclude_deterministic = TRUE`, which left the constants subject to
  variable selection; the element is `exclude_det`. The check now also covers
  the elements of `coef$minnesota`, and rejects `coef$coint_var` for VEC models,
  which do not use it.

* **Documentation for coding assistants, in `inst/agents/`.** An assistant
  working from tutorials writes bvartools code that runs and means something
  else, or that calls functions which no longer exist. `inst/agents/` has an
  `AGENTS.md`, and a skill covering:
  - the workflow and the rules that prevent that: assigning each step back,
    priors given as precisions with no defaults, draws in rows, and
    `vec_to_var()` before forecasting a VEC;
  - the combinations the samplers refuse, and why;
  - complete examples of a VAR, a VEC, a TVP-SV model, a quantile VAR, lag
    order comparison and an HDF5 round trip;
  - references on what `add_priors()` needs for each model type, the layout of
    model objects and their draws, forecasts, impulse responses, sign
    restrictions and spillovers, and in-sample and out-of-sample model
    comparison.

  The installed package carries it at `system.file("agents", package =
  "bvartools")`, matching its version, and the repository is a Claude Code
  plugin marketplace. `tests/testthat/test-agent-docs.R` runs every R example
  in it, and the examples assert the shapes their text states. No function
  changes, so draws are unchanged.

* **Fixes to reading the posterior draws of models whose draws vary by
  period.**

  - `summary()` with `period` reported the error covariance of the *last*
    period for every period, while the coefficients did follow `period`. It
    now reports the covariance of the period asked for.
  - `summary()` and `fevd()` stopped with "subscript out of bounds" on a model
    with time varying parameters and a gamma error term without covariance
    block, because they read its single error precision as one per period. That
    also broke `summary()` on model lists and expanding windows of such models.
    Whether the precision is a path is now read off the stored draws, by
    `summary()`, `plot()` and the draws behind `irf()` and `fevd()` alike,
    where it used to be decided from the specification four different ways.
  - `thin()` left the draws of the state variances of a time varying model
    (`a$sigma`, `psi$sigma`) and of `rho` (`beta$rho`) unthinned while thinning
    the rest, so that rows no longer belonged to the same draw across blocks.
    It now thins every block of draws.
  - `window()` cut the data and left every posterior path at the length of the
    original sample, so that, for instance, the summary of the last period of a
    window was that of a period of the original sample. The paths -- coefficients
    and cointegration vectors of time varying models, error precisions and
    log-volatilities that vary by period, and the pointwise log-likelihood -- are
    now cut to the periods of the window.
  - `selection_criteria()` on an expanding window stopped with "non-numeric
    matrix extent" unless the windows had forecast errors. The in-sample and
    out-of-sample criteria are now each computed when their draws exist, and an
    error is raised only if neither does.

* **New data set `at_macrodata`.** It contains the Austrian sub-model of a
  global VAR model of the 33 countries in the GVAR database of Mohaddes and
  Raissi (2024): quarterly domestic series, their trade weighted foreign
  counterparts and global commodity prices from 1979Q2 to 2023Q3. It was
  produced from data set `gvar2023` of bgvars, and illustrates models with
  weakly exogenous variables. The data set is a list of two time series:
  `domestic` holds the six domestic variables and `foreign` the six foreign and
  three global ones, so they can be passed directly to arguments `data` and
  `exogen` of `create_bvarmodel()` and `create_bvecmodel()`.

* **A time varying cointegration space can be centred on the maximum
  likelihood estimate as well.** For VEC models with `tvp = TRUE`,
  `add_priors(coint = list(rho = 0.999, p_tau_i = "ml", weight = 1))` uses the
  informative marginal prior of Koop et al. (2011, working paper version,
  eq. 12): the transition of the state equation becomes
  `rho (I_r kron P_tau)` with `P_tau = H H' + H_perp T H_perp'`, `H` spanning
  Johansen's estimate of the space, and the state before the sample gets the
  stationary distribution that transition implies. `T` is set so that the
  prior spread of how far beta tilts away from `sp(H)` in each period matches
  the sampling distribution of the estimator, worth `weight` samples. The
  transition is stored as `object$priors$beta$p_tau`.

  `rho` bounds how informative this can be: even `T = 0` leaves a tilt of
  about `sqrt(1 - rho^2)` per period, and `add_priors()` warns when the
  requested weight asks for more.

  A larger weight does not always tighten the posterior. It lowers `T`, and
  `T` is also how much of the tilt persists from one period to the next; near
  `T = 0` each period's tilt is informed by that period's data alone, and the
  posterior can widen again. On `e6` with `rho = 0.999`, `weight = 1`
  (`T` about 0.44) gave a tighter posterior than `weight = 100` (`T = 0`).
  Weights that keep `T` clearly above zero are the useful range.

  A weight so small that `T` is the identity
  in every direction leaves the prior as it was, so such a model's draws are
  those of the noninformative one.

  This needs the updated BayesTS core vendored with it, whose three
  time-varying VEC samplers read `p_tau`; for models without it their draws
  are unchanged.

* **The cointegration space prior can be centred on the maximum likelihood
  estimate.** For VEC models with constant cointegration vectors,
  `add_priors(coint = list(v_i = ..., p_tau_i = "ml", weight = 1))` centres the
  prior of Koop et al. (2010) on the space spanned by Johansen's estimate of
  beta. Its spread is calibrated so that, given the loadings, the prior on how
  far beta tilts away from that space is the estimator's own sampling
  distribution, worth `weight` samples of information. `coint$v_i = "ml"` sets
  the shrinkage of the loadings from their maximum likelihood estimate.
  `coint$p_tau_i` also accepts a full matrix now.

  Such a prior cannot be combined with `coint$v_i = 0`, which is an error: the
  sampler only ever uses the product of `v_i` and `p_tau_i`, so with a zero
  shrinkage any `p_tau_i` gives the uniform prior on the space. And
  `scale_error_correction()` refuses a model that already has a prior with a
  direction, since that direction is stated in the units of the unscaled
  series.

* **A sampler that cannot run now stops.** `add_posterior_coefficients()` used
  to catch the error, print it, and return the *unestimated* model with an
  `error` element added -- an object of the right class, carrying the
  specification and the data of a fitted model and no posterior draws. Nothing
  in this package or in bgvars ever read that flag, so the failure surfaced
  somewhere else entirely: as a missing element several steps later, or as a
  model ranked against its siblings on a log likelihood it did not have. The
  sampler's own message, which says what about the input it could not work
  with, now reaches the caller.

  This changes what a batch does. `add_posterior_coefficients()` on a
  `modellist`, an `expandingwindow` or a bgvars `gvarmodel` is an `lapply` over
  the single-model method, so one unusable specification ends the run instead of
  leaving a hole in the results. If you were relying on a batch completing past
  a model that cannot be estimated, wrap the call in `try()` yourself -- the
  difference being that you then know it happened.

  A user-supplied `posterior_function` is no longer shielded from its own errors
  either.

* **The out-of-sample regressors of a forecast are no longer in SUR form.**
  `add_forecast_input()` stores them as `object$data$forecast$x`, one row per
  forecast period and one column per regressor, where
  `object$data$forecast$z` used to hold the same numbers kroneckered up with an
  identity of order `k` -- `k` times the rows and `k` times the columns.
  `prepare_forecast_input()` returns the matrix under the name `x` for the same
  reason. HDF5 exports write `/data/forecast/x`.

  A forecast applies one `k` by `n` coefficient matrix to one regressor column
  per period, so the wide form held nothing the compact one does not, at `k^2`
  the memory and `k` times the multiplications -- all of the extra ones against
  a structural zero. For the three-variable models in the examples that is
  nothing anyone would notice. It is `k` that decides how much it is worth: a
  174-variable global VAR forecast twelve periods ahead was allocating half a
  gigabyte of regressors, 99.4% of them zeros, and now allocates 17 KB.

  **Draws are unchanged.** A seeded forecast is bit-identical before and after
  for a plain, a structural and a time varying VAR. The vendored core's own
  fixture comparison, which runs larger models under a different BLAS, moves 52
  of 88 fixtures in `/posterior/forecast` and nowhere else, by at most a
  relative 9.4e-16 -- four ulps of reassociated summation.

  **A model fitted with an earlier version still forecasts.** An object or an
  exported file carrying the old `z` is compacted on the way into the sampler,
  exactly -- the kron is subscripted, not averaged -- and a `z` whose
  dimensions are not multiples of `k` is refused rather than turned into
  plausible regressors. Code that reads `object$data$forecast$z` itself has to
  be updated; code that goes through `add_forecast_input()` and
  `add_posterior_forecasts()` does not.

* **The autocorrelation of a time varying cointegration space can be
  estimated.** `add_priors()` takes `coint$rho_min` and `coint$rho_max`, the
  support of a uniform prior on the `rho` of
  `beta_t = rho beta_{t-1} + eta_t`, and with them `VecTvpWishart`,
  `VecTvpGamma` and `VecTvpStochvol` draw `rho` in the Gibbs sampler rather
  than holding it at `coint$rho`. Both ends or neither; `coint$rho` then names
  the value the chain starts at and has to lie between them. The draws come
  back as `object$posterior$beta$rho`, an `mcmc` object of one column. This is
  the block of Koop, Leon-Gonzalez and Strachan (2011) that the package has
  said it was missing since time varying cointegration was added. Nothing
  changes for a model that does not name the two bounds.

  The draw is an exact Gibbs block rather than the Metropolis-within-Gibbs step
  of the paper, and the difference is in the model rather than in the
  algorithm. In the paper the cointegration space at the start of the sample is
  drawn from the state equation's own stationary distribution,
  `N(0, I / (1 - rho^2))`, which makes `rho` appear where no
  conjugacy survives. Here the prior on that state is built once from
  `coint$rho` and does not follow the draw, so it drops out of `rho`'s
  conditional and leaves a normal truncated to the prior's support. The
  shrinkage the loadings carry is fixed the same way. `?add_priors` says so.

* **A time varying cointegration space was centred on the wrong state in the
  first period.** The simulation smoother is handed the prior mean of the state
  the first observation loads on, and it does not put the transition through
  it; `VecTvpWishart`, `VecTvpGamma` and `VecTvpStochvol` were passing
  `beta_0` there unchanged, which is the random walk's answer and is right
  only at `rho = 1`. Below one the smoother was centring `beta_1`
  over `beta_0` while the draw of `beta_0` was centring it over
  `rho beta_0`. **Draws of those three algorithms change** wherever
  `rho` is below one, which is every setting the package suggests; the
  move is about a tenth of a percent at `rho = 0.99` and smaller at the
  0.999 the examples use. No other algorithm is affected. Fixed in the vendored
  BayesTS core, which verified that every other model's draws are unchanged
  digit for digit.

* **Sign restrictions.** `add_sign_restrictions()` identifies the shocks of a
  VAR by the signs of the impulse responses they produce, and `irf()`, `fevd()`
  and `spillover()` use that identification under `type = "sign"`. The
  restrictions are given as a data frame of `impulse`, `response`, `sign` and
  an optional `horizon`, one row per restriction, naming variables rather than
  positions. For every posterior draw the function searches rotations of the
  Choleski factor until one produces the pattern that was asked for; since a
  sign restriction cannot tell a shock from its own negative, a column that
  fails is retried flipped before the rotation is discarded. A draw that no
  admissible rotation is found for within `max_tries` is dropped, and
  `summary()` reports both the restrictions and how much of the posterior
  survived them -- a low rate being a finding about the restrictions rather
  than a technicality. The identification is set valued: what comes back is
  the collection of responses the restrictions admit, not one response per
  draw, so the credible interval spans models rather than estimates. The
  accepted rotations are stored beside the other posterior draws, and are
  carried by `thin()` and by the HDF5 export along with the restriction table.
  Zero restrictions are not supported: they cannot be imposed by rejection and
  need the algorithm of Arias, Rubio-Ramirez and Waggoner (2018). The new
  vignette `sign-restrictions` works through a monetary policy identification
  on `us_macrodata` and says what a set identified credible band does and does
  not mean.

* **`irf()` and `fevd()` return the impact period for `n_ahead = 0`.** A
  horizon of zero is well defined -- the response is `Phi_0 P`, which is `P`,
  and the decomposition rests on `Phi_0 = I` -- but both failed with
  "Mat::cols(): indices out of bounds or incorrectly used" instead. The
  recursions reused the lag order as a count of regressor columns, and on
  impact that count is zero, so the slice of the coefficient matrix asked for
  column `0` through column `-1` on unsigned indices. The count is now separate
  from the lag order and the coefficients are only sliced where the recursion
  actually runs. `irf()` and `fevd()` also validate `n_ahead` the way
  `spillover()` already did, so a negative horizon is rejected in R with a
  message rather than reaching the recursion, and a cumulative impulse response
  over the impact period alone keeps its shape instead of transposing itself.

* **`irf()`, `fevd()` and `spillover()` accept a caller supplied impact
  matrix.** The new type `"custom"` takes the matrix `P` that the forecast
  error responses are post-multiplied by from argument `impact`, either as one
  matrix that identifies every posterior draw the same way or as a list with
  one per draw. Until now the impact matrix could only be named, so every
  identification scheme meant another branch in each of the three recursions;
  it is now an argument, and the recursions need not learn about a scheme to
  carry it. The types that were already there are unchanged, and the matrix
  that reproduces one of them is not always the obvious one -- `irf()`
  normalises the Choleski factor to a unit shock under `"oir"` while `fevd()`
  does not, and the custom path normalises nothing at all. Under `"custom"`
  the two decompositions keep `Sigma` as the forecast error covariance, so
  their shares add up across shocks exactly when `P P'` equals `Sigma`, as it
  does for a rotation of the Choleski factor. This is groundwork for sign
  restrictions, which produce such a rotation per draw.

* **`summary()` no longer marks a covariance that was never estimated as
  significant.** The asterisk says that the credible interval of a value
  excludes zero, and it was placed by comparing the signs of the two bounds.
  Wherever no covariance is estimated -- a gamma prior on the error variances,
  a structural model, a quantile VAR -- both bounds are exactly zero, and zero
  shares its sign with itself, so every off-diagonal entry of the printed
  variance-covariance matrix carried a mark. The test is now the definition:
  the lower bound is above zero, or the upper one below it. Estimated
  coefficients and variances are marked exactly as before, and nothing beyond
  the printed output changes.

* **The structural variance decomposition now weighs each shock by its own
  variance.** `fevd(type = "sir")` used `A_0^-1` as the impulse matrix and
  `A_0^-1 A_0^-1'` as the forecast error covariance, leaving the covariance
  `Sigma` of the structural errors out of both. The decomposition therefore
  treated every structural shock as if it had unit variance, which moved weight
  to whichever shock carried the largest loading in `A_0`: on a US
  inflation/unemployment/interest rate VAR the share of the interest rate
  variance attributed to the unemployment shock came out at 86 per cent against
  the 57 per cent a recursive SVAR of the same model gives. The impulse matrix
  is now `A_0^-1 P`, with `P` the Choleski factor of `Sigma`, and the forecast
  error covariance `A_0^-1 Sigma A_0^-1'`. Decompositions from earlier versions
  are only unaffected where the structural variances happened to be one, so
  published `"sir"` numbers should be regenerated. Impulse responses were never
  affected: `irf(type = "sir")` scales the shock through its own `shock`
  argument, and `shock = "sd"` reproduces `vars::irf()` on the corresponding
  `SVAR` to Monte Carlo error.

* **Reduced-form impulse responses, decompositions and spillovers are refused
  for a structural model.** A structural model stores the contemporaneous block
  separately, so its draws are the structural `A_i` and the covariance of the
  structural errors. The recursion behind `"feir"`, `"oir"` and `"gir"` needs
  the reduced form instead, and reading the structural quantities in its place
  produced responses that belonged to no model at all -- on the same VAR an
  `"oir"` response was out by up to 3.8 against the Choleski response of the
  reduced form. `irf()`, `fevd()` and `spillover()` now stop with an explanation
  rather than return those numbers. Use `"sir"` or `"sgir"` on a structural
  model, or estimate it with `structural = FALSE`.

* **`fevd()` and its `plot()` method can pool the smaller contributions into one
  group.** The new argument `max_groups` caps the number of groups shown: the
  `max_groups - 1` variables with the largest contribution across the whole
  horizon are kept, in the order of the endogenous variables, and the remaining
  ones are added up in a further group `"Other"`. The row sums are untouched, so
  a decomposition that added up to one still does. Passed to `fevd()`, the
  returned decomposition itself has at most `max_groups` columns; passed to
  `plot()`, the full decomposition is kept and only the bars and the legend are
  pooled -- which is the point of the argument, since a model with many
  variables produced a legend that crowded out the bars. The default is `NULL`
  in both, which shows a group per variable as before.

* **Bayesian quantile VARs.** `create_bvarmodel()` takes `error = "ald"` and a
  `quantile`, and estimates a conditional quantile of the endogenous variables
  instead of their conditional mean -- the question behind a growth-at-risk
  exercise, and behind any claim that a relationship differs when things go
  badly. Constant and time varying coefficients both, through the vendored
  `VarNormalAld` and `VarTvpAld` samplers, which arrived with the previous
  refresh and had no way into R until now.

    The model is the asymmetric Laplace one of Kozumi and Kobayashi (2011).
  Minimising the quantile loss at `q` is maximising the likelihood of that
  distribution, and it is a scale mixture of normals, so conditional on the
  latent scales every equation is an ordinary weighted normal regression. That
  is what makes these the stochastic volatility samplers with a different rule
  for where the per-period variance comes from rather than a new kind of model.

    The workflow is the usual one. `add_priors()` takes the inverse gamma prior
  of the scale of the asymmetric Laplace, one per equation, through
  `sigma = list(shape = , rate = )`, and stores it in `priors$u_scale`; note
  that these are used as given, unlike the gamma prior on the error variances,
  which is halved on the way in. `add_initial_values()` starts each scale at the
  mean of the check function of the residuals. `add_posterior_coefficients()`
  and `add_posterior_loglik()` dispatch as they do for every other algorithm,
  and the draws of the scale arrive in `posterior$u_scale`. `summary()` and
  `plot()` report the model as a Quantile-VAR and print the quantile beside the
  lag order.

    **A vector in `quantile` produces a list of models**, one per quantile, in
  the same way a vector of lag orders does -- so a quantile grid is a
  `modellist` that runs through the same functions and can be estimated in
  parallel without the samplers knowing about it.

    **Three things these models do not do**, each because of what a quantile is
  rather than for want of an implementation. They estimate no error covariances,
  since rotating the equations into each other leaves a residual whose quantile
  is not the one that was asked for. They do not forecast, since the `h` step
  ahead quantile is not the quantile of the iterated one step ahead quantiles;
  `add_posterior_forecasts()` says so rather than producing a path that cannot
  be read as one. And the spread of the draws is not a calibrated credible
  interval: the asymmetric Laplace is a working likelihood, the posterior
  locates the quantile, and the adjustment of Yang, Wang and He (2016) is not
  applied. Variable selection is available as `"bvs"`; `"ssvs"` is refused when
  the specification is made.

    A new vignette, *Bayesian Quantile VARs in bvartools*, walks through the
  whole of it, including the check that matters for a model whose failure mode
  is silent: the share of fitted residuals below zero is the quantile. A model
  that had dropped its skew term would estimate the median and pass every other
  test there is.

* **`summary()` and `plot()` label the deterministic terms of every model that
  has them.** Both take those labels from `model$deterministic`, and only
  `create_bvarmodel()` ever recorded it: a model assembled by `bvar()` from the
  draws of a user-written sampler, or converted from a VEC by `vec_to_var()`,
  carried the terms but not their names. The vector of regressor names then came
  back `n` short, and `summary()` stopped with `length of 'dimnames' [2] not
  equal to array extent` instead of labelling anything.

    Both now record the names they already had to hand: the deterministic
  columns `bvar()` isolates from `x`, and the ones `vec_to_var()` collects while
  it rebuilds the regressors of the VAR. Nothing numeric reads these labels, so
  no result changes -- `irf()` and `fevd()` were unaffected throughout, which is
  why this surfaced only in the two methods that print.

    The name helper no longer trusts the field either. A missing or wrong-length
  vector of deterministic names is replaced by `det.1`, `det.2`, ... and of
  exogenous ones by `x1`, `x2`, ..., which is the rule the 'bvecmodel' helper
  already followed -- and the reason no VEC summary ever hit this.

* **Vendored BayesTS core refreshed.** **Draws are unchanged**, for every VAR and
  VEC model this package samples. Upstream's own fingerprint comparison was run
  over the change and reports 78 fixtures recorded before and after, 78
  unchanged and none moved. (Its raw output calls all 78 moved, because the
  recording gained a row -- `/posterior/u_scale/coeffs`, `absent` for every model
  that does not write one. With that row dropped the two recordings are
  identical.)

    What arrived is two models, and nothing else that executes. `VarNormalAld`
  and `VarTvpAld` are Bayesian quantile VARs: they estimate a conditional
  *quantile* rather than a conditional mean, through the normal scale mixture
  representation of the asymmetric Laplace distribution, which makes them the
  stochastic volatility samplers with a different rule for where the per-period
  variance comes from rather than a new kind of model. They are VARs, so they
  belong here rather than to `dfmtools` and are not candidates for the refresh
  script's `skip` list: their samplers, `core/models/ald_support.h` and the
  `core/algorithms/inverse_gaussian.*` draw the latent scales need are all
  compiled into the shared object.

    **Nothing in R can reach them.** There is no `src/VarNormalAld.cpp` or
  `src/VarTvpAld.cpp` binding and no R entry point, so the two sit vendored and
  unreachable -- the position `VecNormalWishart` was in until it got a binding.
  A binding will have to read `VarSpec::quantile`, which `read_spec()` in
  `src/bayests_r_io.h` does not, so the quantile would otherwise arrive as its
  0.5 default and the model would estimate the median while looking like it had
  been asked for something else. It will also have to leave two refusals alone:
  these models take no covariance block and produce no forecast, both rejected
  by `validate()`, because a rotation of the errors is a combination of
  quantiles rather than the quantile of a combination, and an `h` step quantile
  is not the quantile of the iterated one step quantiles.

    Beyond the two samplers the refresh is what they added to the files every
  model shares: their `Initial`, `Input` and `Draws` structs in
  `bayests/inputs.h` and `bayests/results.h`, their `validate()` in
  `core/inputs.cpp`, and `quantile` in `bayests/spec.h`, `VarSpec`'s first
  non-integer member. The only other upstream change to a file this package
  vendors is comment text in `bayests/spec.h` and `core/inputs.cpp` about which
  loading count a FAVAR reads -- a model that is skipped here in any case.

    One thing changed hands rather than arriving. `bayests/arma.h`, the header
  that routes Armadillo through RcppArmadillo, was this package's own and is now
  upstream's, byte for byte, with the destination chosen by a
  `BAYESTS_ARMA_HEADER` define that `src/Makevars` and `src/Makevars.win` set to
  `<RcppArmadillo.h>`. Nothing about the build changes and no file under
  `src/core/` reaches for `<armadillo>` any more, so the rewrite the refresh
  script applies has nothing left to do; it and `src/bayests_rng_guard.cpp`'s
  `static_assert` that Armadillo came out configured for R's RNG both stay, as
  the cheap half of a failure that is otherwise silent.

    `inst/COPYRIGHTS` gained the seven new files. The refresh script compares
  that list against what is vendored and stops until they agree, which is what
  it did here.

* **Writing a model to HDF5 is about a third faster.** `write_to_hdf5()` used to
  ask the file whether a dataset was there before writing it, and then name the
  dataset again for each of its attributes, so a series with two attributes was
  opened three times over. It now writes through the handle that creating the
  dataset hands back. It also closes the handles it opened instead of ending
  with `close_all()`, which finds out what to close by enumerating every open
  object of the file -- an enumeration that cost more per model than writing the
  numbers in it. On a sub-model of 27 models that is 208 ms per model rather
  than 294. What is written is unchanged.

* **`write_to_hdf5()` records the chain of the `psi` draws correctly.** The
  `start`, `end` and `thin` of the draws in `posterior/psi` were attached to the
  datasets of `posterior/a` rather than to their own, so `psi` came back without
  the parameters of its chain, and those of `a` were overwritten with `psi`'s.
  The blocks that write each set of draws into a group of its own are now one
  loop over the names, which is what let two of them disagree.


* **A model can be addressed by group, so one HDF5 file can hold several.**
  `write_to_hdf5()` and `read_model_from_hdf5()` take a `group` argument, and
  the new `list_models_in_hdf5()` reports which groups of a file hold a model.
  Without a `group` nothing changes: the model is the whole file, and an
  existing file is still refused.

    A model is a group with a `model` subgroup carrying an `algorithm`
  attribute, and the search stops at one rather than descending into its
  `data`, `priors` and `posterior`. The spelling of a group name -- leading
  slash, no trailing slash, `""` for the root -- is the one the BayesTS command
  line uses, so a group named here and a `--group` or `--all-groups` passed
  there mean the same thing. Verified against it: R and `bayests` list the same
  three models of a file that also holds a non-model group, and neither picks
  up the latter.

    What must not already exist is the model rather than the file, which is what
  lets a second model be written beside the first. A failure part way through
  undoes only what the call created: a file it made is removed, a group it added
  to an existing file is unlinked on its own, so the models beside it survive.

* **`read_models_from_folder()` returns a flat, named list.** It used to mirror
  the directory tree into a nesting and give the result no names at all, so a
  folder of sub-models came back as `[[1]]`, `[[2]]`, `[[3]]` with no way to
  tell which was which except by re-deriving it from the order `list.dirs()`
  happened to return. Each model is now named for the file it came from,
  relative to the folder, with its group appended where a file holds more than
  one. Files holding several models contribute all of them.

    The class of the collection is read off the models instead of being guessed
  from a path. An expanding window was recognised by looking for `"ExpWind"` in
  the first file's *full* path, so the class of the result depended on the names
  of every directory above it as well. `write_to_hdf5()` now records the
  collection in the models themselves; the old name is still honoured, so
  folders written before this are read the same way as before.

* **`write_to_hdf5()` no longer swallows a failed write.** The bodies of the
  'bvarmodel' and 'bvecmodel' methods were wrapped in `try()` whose result was
  discarded, so any failure part way through -- an unwritable path, a full disk,
  an element of `object` HDF5 cannot store -- was silently absorbed. The
  function returned as though it had worked and left a half-written file on
  disk. Worse, the "File already exists" guard then turned the *next*, correct
  attempt at the same path into a second and misleading error.

    The error now reaches the caller. On the way out of a failure the HDF5
  handle is closed, so the file is not left locked for the rest of the session,
  and the incomplete file is removed, so a retry meets the real problem rather
  than the leftovers. Only a file the call created is removed; an existing one
  is still refused before the handle is opened, and is left untouched.

    A successful write now returns the path invisibly instead of the open
  `H5File` object it used to leak.

    **Draws are unchanged.** No sampler is involved: this is the file writer
  only. The whole test suite passes, and `test-hdf5.R` gained six tests covering
  the failure path, which the existing round trips never reached.

* **Vendored BayesTS core refreshed.** **Draws are unchanged**, for every VAR and
  VEC model this package samples. Upstream's own fingerprint comparison was run
  over the change and reports 76 fixtures unchanged and none moved, twice: once
  over the shared algorithm it touched and once over the sampler added on top of
  it.

    Nothing that arrived is reachable from here. Upstream added
    `FavarNormalWishart`, a factor augmented VAR, which is a factor model and so
    belongs to `dfmtools`; it is in the refresh script's `skip` list beside the
    four dynamic factor models, along with `core/models/favar_support.h`. That
    file includes the already-skipped `core/models/dfm_support.h`, so an
    unlisted copy would fail to compile rather than sit unreachable -- the same
    mechanism the previous refresh describes, and the reason the list has to be
    extended by hand.

    What does arrive is what the new model added to the files every model shares:
    `FavarNormalWishartInitial`, `FavarNormalWishartInput` and
    `FavarNormalWishartDraws` in `bayests/inputs.h` and `bayests/results.h`, its
    `validate()` in `core/inputs.cpp`, `n_obs_factors` with `n_state()`,
    `n_favar_lambda()` and `n_favar_a()` in `bayests/spec.h`, and
    `TrainData::f_obs` in `bayests/data.h`. The same arrangement the DFMs have
    had since they moved to `skip`: their type surface compiles in, their
    samplers do not.

    One thing beyond that, in a file this package does vendor.
    `core/algorithms/chan_jeliazkov_2009.cpp` gains a second entry point,
    `chan_jeliazkov_2009_conditional`, which holds the trailing elements of every
    state column at observed values instead of drawing them -- what a state
    vector that is part data needs. Nothing here calls it. Its arrival split the
    existing function into an assembly, a conditioning step and a draw; upstream
    verified that split moved no fingerprint, and the three time-varying models
    here that reach the band sampler are unaffected.

* **`create_external_forecast()`** puts forecasts that were produced elsewhere
  into a horse race against the models of this package. It takes any data frame
  of point forecasts in long format -- a period, a publication date, a variable
  and a value, under whatever column names the source uses -- and returns an
  object that mimics the windows of `use_expanding_window()`. From there the
  usual path applies: `combine_models()` puts it beside the models,
  `add_forecast_errors()` scores it against the same test data, and
  `selection_criteria()`, `print(relative = )`, `plot()` and
  `plot_forecast_errors_by_period()` treat it as one more competitor. The
  functions that add priors, initial values or posterior draws are without
  effect for it, so a list that mixes models and external forecasts can be run
  through the workflow in one go.

    An external forecast has a publication date where a model has a training
  sample, and the two are not the same thing: the forecaster did not know the
  observations that had not been released yet. Each publication is therefore
  matched to the last training sample that ends no later than `data_lag` periods
  before it, which defaults to one period. With annual data that puts a forecast
  published during 2021 against a model estimated through 2020, so the forecast
  for 2021 is the one-step ahead forecast of both. Where several publications
  fall to the same training sample -- two forecast rounds in one year, say --
  only one is used, the latest by default, because a forecaster must not enter
  the same comparison twice.

    External forecasts are point forecasts, so each publication contributes a
  single value rather than a posterior: the credible bands of their out-of-sample
  statistics are degenerate and they have no in-sample criteria. Three functions
  that previously assumed every model carries both kinds of criteria --
  `print.selcritlist()`, `plot.selcritlist()` and `choose_best_model()` -- now
  leave a model's row empty instead of failing when a criterion is missing from
  it.

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

* **Fixed: `irf()` had both defects the `fevd()` fix above describes.** **Draws
  are unchanged; impulse responses of time-varying, stochastic volatility and
  structural models change, and the old numbers were wrong.** `.collect_draws()`
  was introduced for `fevd()` and `spillover()` but `irf.bvarmodel()` kept its
  own copy of the code it replaced, so it carried the same `nrow(train$y) / k`
  sample length -- reading a truncated row for the default `period` and refusing
  any `period` past `tt / k` -- and, additionally, had the line that folds
  `A_0^{-1}` into the coefficients commented out, so `type = "sir"` and
  `type = "sgir"` built the recursion from the structural `A_i` instead of the
  reduced-form `A_0^{-1} A_i` that `fevd()` uses. Constant-coefficient,
  non-structural models are unaffected: they read neither `tt` nor `A_0`.

    `irf.bvarmodel()` now calls `.collect_draws()` like the other two, and only
    attaches the shock size itself. The helper gained a `need_Sigma` argument so
    that a forecast error or structural response does not pay for one covariance
    inversion per draw to obtain a matrix it never reads.

* **Fixed: `create_bvecmodel(seasonal = ...)` failed on data of frequency one.**
  The branch that builds the dummies is skipped at a frequency of one, with a
  warning, but the `cbind()` calls that add them to the error correction term or
  to the unrestricted regressors sat outside it and reached for a `seas` that was
  never created, so the call stopped with `object 'seas' not found` instead of
  proceeding without seasonal terms. They now sit inside the branch, as they
  already did in `create_bvarmodel()`.

* **Fixed: the Minnesota prior computed a least squares estimate it did not
  need, and failed on short samples because of it.** `minnesota_prior()`
  obtained the residual covariance of the full VAR form before choosing between
  `sigma = "AR"` and `sigma = "VAR"`, though only the latter ever reads it. That
  estimate needs more observations than the model has regressors per equation;
  the default does not, since it regresses each variable on its own lags and the
  deterministic terms alone, which is a far smaller system. A model with more
  regressors per equation than training observations therefore stopped at a
  singular matrix from `solve()` under the default, for a quantity that was then
  discarded -- five variables at four lags is twenty-one coefficients, which a
  twenty-five quarter sample cannot carry once the lags are taken off, while
  each of the univariate regressions on that same sample has five regressors
  against twenty observations and is perfectly well conditioned. The estimate is
  now obtained only where `sigma = "VAR"` asks for it. Priors are unchanged
  wherever the old code ran, on both paths and for VAR and VEC models alike.

    Both paths now also check that the sample can support the regression they
    are about to run, and name the two counts. The VAR path needed that for more
    than the message: `tt - nrow(x)` is the denominator of the covariance, so a
    sample exactly as wide as the regressor set divided by zero and a narrower
    one by a negative number, which returned an infinite or sign-flipped prior
    variance wherever the inverse happened to succeed rather than failing at
    all. `ssvs_prior()` gained the same check under `semiautomatic`, where the
    least squares standard errors are genuinely used and the requirement
    therefore stands; the fixed `tau` values remain available on a sample that
    cannot meet it.

* **Fixed: `ssvs_prior()` did not work on a `bvecmodel` at all.** The function
  read `object$data$y`, `$w`, `$x` and `$z`, the layout that preceded the move of
  the estimation sample under `object$data$train`, so every one of them was
  `NULL` and the call stopped at `t(NULL)` with "argument is not a matrix". Its
  own example failed. `ssvs_prior.bvarmodel()` was already on the current layout.

    Reaching the rest of the function exposed a second defect: under
    `semiautomatic`, `tau1` was assembled by appending to `tau0` rather than to
    itself, so it came out with the alpha block, the whole of `tau0` and then its
    own values -- longer than the coefficient vector it belongs to, and wrong
    where it overlapped. Both vectors now have `ncol(z)` elements for every rank
    and lag order, and `tau1` exceeds `tau0` element by element as the
    semiautomatic approach intends.

* **Fixed: four more places where the VEC code read the pre-`train` data
  layout.** Unlike `ssvs_prior()` above, all four sat in an `is.null()` guard, so
  they did not fail -- they took the wrong branch in silence.
  `inclusion_prior.bvecmodel()` tested `object$data$z` while its body already
  read `object$data$train$z`, so the body never ran and the function returned
  `NULL`, which surfaced downstream as "a prior inclusion probabilities must have
  n elements, got 0". `.check_bvecpost_input()` skipped every coefficient-side
  check for VEC models. In `add_priors.bvecmodel()`, `coef$const` was never
  applied and the least squares covariance the Minnesota prior needs for its
  analytical solution was never stored.

    Bringing each branch to life exposed what had drifted while it was
    unreachable. `inclusion_prior.bvecmodel()` defined `n_gamma` and `n_upsilon`
    inside its Minnesota-like branch but located the deterministic block from
    them afterwards, so the default `minnesota_like = FALSE` would have stopped
    at "object 'n_gamma' not found"; they are now defined once for both, as in
    `inclusion_prior.bvarmodel()`. The `coef$const` block added its column offset
    `r` twice, once to the position and once again when indexing, and had lost
    the `is.numeric()` guard around its numeric branch, so `coef$const = "mean"`
    would have written the string into the prior mean and coerced the whole
    matrix to character. `.check_bvecpost_input()` still asked for `v_i`,
    `a_v_i` and `psi_v_i`, the names these elements carried before they became
    `v_inv`, `a_sigma_inv` and `psi_sigma_inv`; it now matches
    `.check_bvarpost_input()` element for element.

* **SSVS and BVS now work on VEC models.** `add_initial_values.bvecmodel()` had
  no counterpart to the `a_lambda` and `psi_lambda` block of
  `add_initial_values.bvarmodel()`, so a VEC model with `varsel` set stopped at
  "a initial inclusion indicators must have n elements, got 0". It now writes
  both, one indicator per element of the coefficient vector rather than one per
  selected position: the sampler masks the regressors with `diag(lambda)` as a
  whole, so a position the sweep never visits keeps whatever it starts with for
  the whole run. Every position therefore starts at one, which leaves the
  unselected coefficients in the model -- in particular the `k * r` loadings at
  the front of `a`, which a zero would have switched off permanently.

    `inclusion_prior.bvecmodel()` had to be corrected alongside it. It built the
    list of excluded positions correctly, loadings first, but applied it only
    inside `if (n_c_unres > 0 & exclude_deterministics)`, and `varsel$exclude_det`
    defaults to `FALSE`. So in the default case the loadings stayed in `include`
    and the sampler refused them: "variable selection cannot be applied to the
    k * r loading coefficients at the front of a VEC's a". The exclusion is now
    applied whether or not the deterministics are also dropped.

    Verified for SSVS and BVS at ranks one and two, with and without
    unrestricted deterministic terms, and for `error = "gamma+covar"` with
    `varsel$covar = TRUE`, which exercises `psi_lambda`. Across draws the
    loadings stay switched on and the selected positions move.

    Unrelated to variable selection, `add_initial_values.bvecmodel()` read
    `object$priors$a$v_i` in its `method = "prior"` branch, one more of the
    pre-rename names; it is now `v_inv`, as it already was thirty lines below for
    `psi`.

* **Variable selection alongside an error covariance block is refused where the
  sampler cannot do it, instead of failing inside the sampler.** A
  constant-coefficient model with `error = "gamma+covar"` or `"sv+covar"` and
  `varsel` set, but without `varsel$covar`, used to reach the sampler and stop
  there on "psi prior inclusion probabilities must have n elements, got 0" -- a
  prior it was never given, named after an argument the caller had not set.

    The samplers for those models read one selection scheme for the whole model,
    so a covariance block they are given is selected along with the
    coefficients; only the time-varying ones take the covariance block's scheme
    separately, which is why `varsel$covar = FALSE` is a real choice there and
    not here. That asymmetry lives in the vendored core -- the constant
    coefficient input types have no per-block field to set -- so `add_priors()`
    now says so up front, and names the three ways out: select the covariances
    too, drop them, or use a time-varying model. `varsel$covar = TRUE` and every
    combination without a covariance block are unaffected.

* **Fixed: maximum likelihood initial values of a rank zero VEC model were
  fitted against the wrong response.** **Draws from such a model change, and the
  old starting values were arbitrary.** `add_initial_values.bvecmodel()`
  transposed `y` inside its `if (r > 0)` block, but the least squares fit below
  stacks it with `matrix(y)` whether or not there is a cointegration term to
  estimate first. With no such term the series went in stacked variable by
  variable where the SUR regressors want it period by period, so the fit was
  against a permuted response. Nothing failed -- the chain simply started
  somewhere unrelated to the data. This is the default `method` and
  `create_bvecmodel()` generates `r = 0` among its default ranks, so it was
  easy to reach. The transpose now happens for every rank, and the initial
  values reproduce the least squares fit their regressors imply at ranks zero,
  one and two.

    The fallback residual carried the same confusion: `matrix(y, k)` on the wide
    series fills each column with consecutive observations of one variable
    rather than one period across variables. It is now transposed. That value
    survives only where the sample is too short for the fit above, the path that
    warns and sets the coefficients to zero, and it is what the error
    covariance's own initial value is then built from.

* **Fixed: `add_initial_values(method = "prior")` on a VEC model.** The branch
  transposed the error correction regressors but not the endogenous variables,
  and then formed the residual against the SUR regressors, which need the series
  stacked variable-within-period. The subtraction was non-conformable and the
  call stopped at "non-conformable arrays". Both are now transposed, as the
  maximum likelihood branch above already did, and the residual is only formed
  where there are regressors to form it from. Verified for ranks one and two,
  with and without unrestricted deterministic terms, and for `gamma+covar` and
  time-varying specifications; each now yields initial values the sampler
  accepts.

    `method = "prior"` on a VAR model turned out not to be broken. What looked
    like a defect was an improper prior: drawing Sigma from a Wishart with fewer
    degrees of freedom than endogenous variables is not a draw from anything,
    and `stats::rWishart()` refuses it with a message that names neither the
    argument nor the function that set it. Since `sigma$df` below `k` is a
    perfectly good prior to *estimate* with -- the posterior adds one degree of
    freedom per observation, which is why the fixtures use `df = 1` -- the
    restriction belongs to `method`, not to the prior, and is now reported that
    way, pointing at `sigma$df` and at `method = "maxlik"`. `maxlik` is
    unaffected.

* **Fixed: three examples that no longer matched the API.** `predict.bvarmodel()`
  and `plot.bvarprd()` called `predict()` straight after
  `add_posterior_coefficients()`, from before forecasting moved into
  `add_forecast_input()` and `add_posterior_forecasts()`; they now run those two
  first. `post_coint_kls()` read `temp$data$y` and passed `w = ect` where the
  example only ever defined `w`; it now follows `post_coint_kls_sur()`, which had
  already been migrated, and sizes the priors from the data rather than by a
  literal. With these and the missing comma above, all 88 documented examples run.

* **`NAMESPACE` is no longer generated from a blanket pattern.** The package
  carried `exportPattern("^[[:alpha:]]+")` from its Rcpp skeleton alongside an
  otherwise explicit export list. It was exporting three internal plotting
  helpers marked `@noRd` -- `draw_error_bars()`, `forecast_plot_series()` and
  `plot_selcrit_forecast_errors()` -- which is what `R CMD check` reported as
  undocumented code objects; those are now internal. In the other direction,
  fifteen documented compiled functions (`post_normal()`, `post_normal_sur()`,
  `post_bvs()`, `post_coint_kls()`, `post_coint_kls_sur()`,
  `post_gamma_measurement_variance()`, `post_gamma_state_variance()`, `ssvs()`,
  `coint_prepare_sur_data()`, `coint_kls2010_reparameterise_two()`,
  `covar_prepare_data()`, `covar_vector_to_matrix()`,
  `generate_lower_block_diagonal()`, `loglik_normal()` and
  `sur_const_to_tvp()`) reached the namespace only through that pattern and now
  carry `@export` in their own source, as `kalman_durbin_koopman_2002()` and the
  two stochastic volatility functions already did. The set of exported objects is
  otherwise unchanged.

* **Fixed: a missing comma in the examples of `generate_artificial_var()`.** It
  left `bvartools-Ex.R` unparseable, so `R CMD check` ran no examples at all.

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

# bvartools 0.3.0

This is a transition release. It is the functionality of the previous CRAN
version with the fixes listed below, and it is the last release before 1.0.0
reorganises the package around a different set of functions. Nothing here stops
working.

* **Functions that bvartools 1.0.0 does not have any more announce themselves.**
  The first time in a session that such a function is used it emits a message
  naming its successor, or saying that it has no replacement. The message is
  shown once per function per session and can be switched off with
  `options(bvartools.transition.messages = FALSE)`. It is a message rather than
  a warning, so it cannot become an error under `options(warn = 2)`.

    Renamed: `gen_var` to `create_bvarmodel`, `gen_vec` to `create_bvecmodel`,
  `bvec_to_bvar` to `vec_to_var`, `kalman_dk` to `kalman_durbin_koopman_2002`,
  `stochvol_ksc1998` to `stochvol_ksc_1998`, `stochvol_ocsn2007` to
  `stochvol_ocsn_2007`, `stoch_vol` moved into `stochvol_ksc_1998`, and `bvs` to
  `post_bvs`. Replaced by a different workflow: `draw_posterior`, `bvarpost`
  and `bvecpost`, which become `add_posterior_coefficients` alongside
  `add_posterior_forecasts` and `add_posterior_loglik`. Removed with no
  successor: `post_normal_covar_const`, `post_normal_covar_tvp`, and the whole
  dynamic factor model branch -- `dfm`, `dfmpost`, `gen_dfm`, their methods, and
  the data set `bem_dfmdata`, which cannot announce itself.

    Methods on the renamed classes `bvar`, `bvec` and `bvarlist` are silent,
  because the generic that dispatches them is unchanged and only the class is
  renamed, to `bvarmodel`, `bvecmodel` and `modellist`. So are the functions
  that keep their name in 1.0.0 but take the reorganised model object:
  `add_priors`, `bvar`, `bvec`, `irf`, `fevd`, `inclusion_prior`,
  `minnesota_prior` and `ssvs_prior`. The new vignette, *Moving from bvartools
  0.3.0 to 1.0.0*, lists all of it.

* **Fixed: `stochvol_ksc1998` and `stochvol_ocsn2007` failed on an observation
  far out in the tails of every mixture component.** Both sampled the mixture
  indicator from weights formed as densities and normalised by their sum. Where
  the log of the squared observation lies far enough below the log-volatility --
  about 112 for the seven components of Kim, Shephard and Chib, about 158 for
  the ten of Omori, Chib, Shephard and Nakajima -- every density underflows to
  zero, the row sums to zero, the weights become `NaN`, and the sampled
  indicator runs one past the last component, ending the call with
  `Mat::elem(): index out of bounds`. The weights are now formed in logs and
  shifted by their row maximum before they are exponentiated, and the indicator
  is clamped to the components that exist.

    This is algebraically the same calculation, and draws are unchanged where
  they were being produced at all: verified bit for bit against the previous
  implementation over the `us_macrodata` series from a fixed seed. `stoch_vol`
  is a wrapper for the first of the two and inherits the fix.

* **Fixed: neither function checked the size of `sigma`, `h_init` or
  `constant`.** They were indexed on trust, so a vector of the wrong length was
  reported as `Mat::elem(): index out of bounds` instead of as a statement about
  the argument. Each is now checked against the number of columns of `y`.

* **Fixed: `irf` and `fevd` produced reduced form quantities from a structural
  model.** A structural model keeps its contemporaneous block separately, so its
  coefficient draws are the structural `A_i` and its covariance draws the
  covariance of the structural errors. The forecast error, orthogonalised and
  generalised recursions want the reduced form. Given the structural quantities
  they returned numbers that belong to no model at all, and the same numbers for
  all three types, since the structural error covariance makes the
  orthogonalisation degenerate. `irf` with `type` of `"feir"`, `"oir"` or
  `"gir"`, and `fevd` with `"oir"` or `"gir"`, now stop on a structural model
  and say why. The structural types `"sir"` and `"sgir"` are unchanged, as is
  every reduced form model.

* **Fixed: the structural variance decomposition ignored the variances of the
  structural shocks.** `fevd` with `type = "sir"` used `A_0^-1` as the impulse
  matrix and `A_0^-1 A_0^-1'` as the forecast error covariance, leaving the
  covariance of the structural errors out of both, so every structural shock was
  decomposed as if it had unit variance. Weight moved to whichever shock loads
  most heavily in `A_0`, and since the shares are normalised they still summed
  to one, so the output gave nothing away. On a three variable example whose
  shock variances stand at 4, 1 and 0.25 the reported shares were 0.21, 0.13 and
  0.65 where they should be 0.74, 0.12 and 0.14. The impulse matrix is now
  `A_0^-1 chol(Sigma)'` and the forecast error covariance `A_0^-1 Sigma A_0^-1'`.
  `"sgir"` already carried `Sigma` and is unchanged, as is every reduced form
  type.

* **Fixed: `gen_vec` stopped on seasonal terms for data of frequency one.** It
  warned that no seasonal dummies are generated and then added them anyway,
  failing with `object 'seas' not found` because the dummies had never been
  built. It now does what the warning says. `gen_var` was never affected.

* Added a test suite, which the package did not have before. It covers the
  announcements and carries one regression test per fix above.


* Carried over from the development version that was never released as 0.2.5.
  All six survive into 1.0.0 and so announce nothing.

    * Added function `covar_vector_to_matrix`.
    * Added function `sur_const_to_tvp`.
    * Updated `Rcpp` dependency in DESCRIPTION file to version 1.0.12.
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
