# Convergence Diagnostics of Several Chains

Compares the chains of a model estimated with `chains` above one in
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
and reports, for every parameter, the potential scale reduction factor
and two effective sample sizes.

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
`block`, `parameter`, `rhat`, `ess_bulk` and `ess_tail`.

## Details

A single chain can look converged and not be: its autocorrelations and
its effective sample size describe how well it explores the region it is
in, not whether that region is the posterior. Chains started from the
same values but drawing different random numbers settle in different
places if the posterior has several modes or the sampler has not left
the neighbourhood of its start, and comparing them is the check a single
chain cannot provide.

The statistic is the rank-normalised split \\\hat{R}\\ of Vehtari et al.
(2021). Every chain is cut into two halves, the pooled draws are
replaced by their ranks and put back on a normal scale, and the variance
between the means of the halves is compared with the variance within
them, \$\$\hat{R} = \sqrt{\frac{\frac{n - 1}{n} W + \frac{1}{n}
B}{W}},\$\$ where \\n\\ is the length of a half, \\W\\ the average
variance within the halves and \\B\\ \\n\\ times the variance of their
means. Splitting the chains catches a chain that is still drifting.

**The rank normalisation is what makes the number mean the same thing
whatever scale the parameter is on.** The plain version of this
statistic is built on variances, so it is undefined for a posterior
heavy-tailed enough to have none and it moves when the same draws are
transformed monotonically – which is awkward for a package whose
parameters are freely rescaled and whose variance blocks are reported
both as `omega` and as its square. What is reported is the larger of the
rank-normalised \\\hat{R}\\ and the one computed on the draws folded
around their median, so that chains agreeing about the centre but not
about the spread are caught as well.

Values close to one indicate that the chains describe the same
distribution; above 1.01, Vehtari et al. (2021) recommend running the
chains longer or reconsidering the model. A parameter that does not vary
in any chain, such as a coefficient that variable selection excluded
throughout or an equation that `iid` left without coefficients, has no
\\\hat{R}\\ and is reported as `NA`.

The signed standard deviations `omega` of a block estimated under the
non-centred prior `omega_v` are symmetric around zero by construction –
the sampler switches their sign at random – so their \\\hat{R}\\ is
close to one whether or not the chains agree. Their squares, `sigma`,
which are in the posterior beside them, are the ones to read.

**Two effective sample sizes are reported, and for this package the
second is usually the one that matters.** `ess_bulk` is computed on the
rank-normalised draws and says how much independent information the
sample carries about the centre of the posterior. `ess_tail` is the
smaller of the effective sample sizes at the 5th and 95th percentiles
and says the same about the extremes. Almost nothing this package
reports is a point:
[`irf`](https://franzmohr.github.io/bvartools/reference/irf.md),
[`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md) and
the forecasts all come back as quantiles, and it is `ess_tail` that says
whether those quantiles have settled. A sample can carry a comfortable
`ess_bulk` and still have a band that moves from one run to the next.

Neither is a convergence diagnostic: chains stuck in different modes can
each have a large one.

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
#> [1] 1.009867
```
