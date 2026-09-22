# Convergence Diagnostics of Several Chains

Compares the chains of a model estimated with `chains` above one in
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
and reports, for every parameter, the split potential scale reduction
factor and the effective sample size.

## Usage

``` r
chain_diagnostics(object, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel' or 'bvecmodel', whose posterior was
  simulated with `chains` of at least two.

- ...:

  further arguments passed to or from other methods.

## Value

A data frame with one row per parameter of every block of posterior
draws, other than the forecasts and the log-likelihood, and the columns
`block`, `parameter`, `rhat` and `ess`.

## Details

A single chain can look converged and not be: its autocorrelations and
its effective sample size describe how well it explores the region it is
in, not whether that region is the posterior. Chains started from the
same values but drawing different random numbers settle in different
places if the posterior has several modes or the sampler has not left
the neighbourhood of its start, and comparing them is the check a single
chain cannot provide.

The statistic is the split \\\hat{R}\\ of Gelman et al. (2013, section
11.4): every chain is cut into two halves, and the variance between the
means of the halves is compared with the variance within them,
\$\$\hat{R} = \sqrt{\frac{\frac{n - 1}{n} W + \frac{1}{n} B}{W}},\$\$
where \\n\\ is the length of a half, \\W\\ the average variance within
the halves and \\B\\ \\n\\ times the variance of their means. Splitting
the chains also catches a chain that is still drifting. Values close to
one indicate that the chains describe the same distribution; above 1.01,
Vehtari et al. (2021) recommend running the chains longer or
reconsidering the model. A parameter that does not vary in any chain,
such as a coefficient that variable selection excluded throughout, has
no \\\hat{R}\\ and is reported as `NA`.

The signed standard deviations `omega` of a block estimated under the
non-centred prior `omega_v` are symmetric around zero by construction –
the sampler switches their sign at random – so their \\\hat{R}\\ is
close to one whether or not the chains agree. Their squares, `sigma`,
which are in the posterior beside them, are the ones to read.

The effective sample size is the sum over the chains of
[`effectiveSize`](https://rdrr.io/pkg/coda/man/effectiveSize.html) of
each chain. It is not a convergence diagnostic: chains stuck in
different modes can each have a large one.

The chains start from the same initial values, those of
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md),
and differ in their random numbers. Dispersed starting values would make
the comparison more demanding.

## References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian data analysis* (3rd ed.). Boca Raton: CRC
Press.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved
\\\hat{R}\\ for assessing convergence of MCMC. *Bayesian Analysis,
16*(2), 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_coefficients.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_files()`](https://franzmohr.github.io/bvartools/reference/bayests_files.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r
data("e1")
e1 <- diff(log(e1)) * 100

model <- create_bvarmodel(e1, p = 1, deterministic = "const",
                          iterations = 200, burnin = 100)
model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = "k", scale = 1))
model <- add_initial_values(model)
model <- add_posterior_coefficients(model, chains = 2)

diag <- chain_diagnostics(model)
max(diag$rhat, na.rm = TRUE)
#> [1] 1.009975
```
