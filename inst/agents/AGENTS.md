# AGENTS.md — using bvartools

Guidance for coding agents writing R code with
[bvartools](https://github.com/franzmohr/bvartools), which estimates Bayesian VAR
and VEC models. The full skill is `skills/bvartools/SKILL.md`, with tested
references in `skills/bvartools/references/`: `recipes.md` for complete
examples, `priors.md` for what `add_priors()` needs, `objects.md` for the layout
of models and draws, `analysis.md` for forecasts and impulse responses, and
`comparison.md` for model selection. This file is the part worth keeping in
context at all times.

1. **Assign every step back**: `model <- add_priors(model, ...)`. Each function
   returns the model with something added.
2. **The order**: `create_bvarmodel()` or `create_bvecmodel()`, then
   `add_priors()`, `add_initial_values()` and `add_posterior_coefficients()`.
   Forecasting adds `add_forecast_input()` and `add_posterior_forecasts()`, then
   `predict()`. Model comparison adds `add_posterior_loglik()`, then
   `selection_criteria()`.
3. **`add_priors()` has no defaults.** `coef` and `sigma` are required, `coint`
   for a VEC, and prior variances are precisions. Which `sigma` elements are
   needed depends on `error`: `df` and `scale` for the default Wishart prior,
   `shape` and `rate` for gamma, six elements for stochastic volatility. A
   missing or misspelt element stops with an error. Check `model$priors`
   afterwards.
4. **Posterior draws are rows**: `model$posterior$<block>$coeffs` is a
   `coda::mcmc` matrix, draws × parameters. A posterior mean is `colMeans()`.
5. **A VEC is analysed in levels**: `vec_to_var()` before `irf()`, `fevd()` or
   `spillover()`. `add_forecast_input()`, `add_posterior_forecasts()` and
   `predict()` take the `bvecmodel` itself and forecast in levels.
6. **Vector arguments** (`p = 1:3`) create a `'modellist'`, and every step maps
   over it.
7. **Refusals are statistical**: no structural model with a Wishart prior, no
   SSVS with stochastic volatility or time-varying parameters, and no forecasts
   or covariances for a quantile VAR (`error = "ald"`).
8. **`gen_var()` and `gen_vec()` are gone.** Tutorials that use them, or
   `$data$Y` and `$data$SUR`, predate the current API.
9. **Keep `set.seed()`**: it reaches the C++ samplers, so a seeded run
   reproduces. If two seeds give different summaries, the chain is too short:
   create the model with `thin` to run it longer without storing more draws.
   `add_posterior_coefficients(model, chains = 4)` runs that comparison for you,
   and `chain_diagnostics()` reports its split R-hat; a single chain cannot show
   that it is stuck in one mode.
10. **Read the method's help page**, e.g. `?add_priors.bvarmodel` against
    `?add_priors.bvecmodel`, rather than guessing a list element.

The installed package carries these files at
`system.file("agents", package = "bvartools")`, matching its version.
Dynamic factor models are in `dfmtools`, which has its own guide.
