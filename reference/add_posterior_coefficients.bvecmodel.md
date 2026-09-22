# Posterior Simulation of Model Coefficients

Forwards model input to posterior simulation functions for vector error
correction models.

## Usage

``` r
# S3 method for class 'bvecmodel'
add_posterior_coefficients(
  object,
  posterior_function = NULL,
  chains = NULL,
  ...
)
```

## Arguments

- object:

  an object of class 'bvecmodel', usually, a result of a call to
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
  in combination with
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
  and
  [`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md).

- posterior_function:

  the function to be applied to the model in argument `object`. If
  `NULL` (default), internal functions are used.

- chains:

  the number of chains to simulate. If `NULL` (default), the value in
  `object$model$chains` is used, and one chain when there is none. Each
  chain is the same simulation with a seed of its own – the first with
  the seed of the model, so that one chain draws what it always has –
  and the chains are pooled, one after the other, in the draws of
  `posterior`, so that every later step uses all of them. Their number,
  when above one, is stored in `object$model$chains`, and
  [`chain_diagnostics`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md)
  compares them. A `posterior_function` is called once per chain. Not
  available for discounted models, whose posterior is not a chain.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the element `posterior` added, whose
elements hold the draws after burn-in as
[`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) objects with one row
per draw and one column per parameter, in element `coeffs`: `beta`, the
cointegration coefficients, \\K\_\beta r\\ columns or \\T K\_\beta r\\
for time varying cointegration, `a`, the loadings and the remaining
coefficients, and the draws of the error term as described in
[`add_posterior_coefficients.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md).

## Details

The function implements commonly used posterior simulation algorithms
for Bayesian VAR models with both constant and time varying parameters
(TVP) as well as stochastic volatility. It can produce posterior draws
for standard BVAR models with independent normal-Wishart priors, which
can be augmented by stochastic search variable selection (SSVS) as
proposed by Geroge et al. (2008) or Bayesian variable selection (BVS) as
proposed in Korobilis (2013). Both SSVS or BVS can also be applied to
the covariances of the error term.

The implementation follows the descriptions in Chan et al. (2019),
George et al. (2008) and Korobilis (2013). For all approaches the SUR
form of a VAR model is used to obtain posterior draws. The algorithm is
implemented in C++ to reduce calculation time.

The function also supports structural BVEC models, where the structural
coefficients are estimated from contemporary endogenous variables, which
corresponds to the so-called (A-model). Currently, only specifications
are supported, where the structural matrix contains ones on its diagonal
and all lower triangular elements are freely estimated. Since posterior
draws are obtained based on the SUR form of the VEC model, the
structural coefficients are drawn jointly with the other coefficients.

The internal samplers draw with the seed in `object$model$seed`, which
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
sets and
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
replaces. R's random number generator is set to that seed, with R's
default kinds, for the simulation and put back as it was afterwards. A
call of [`set.seed()`](https://rdrr.io/r/base/Random.html) between
[`add_initial_values()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
and this function therefore does not change the draws. A model without a
seed draws from R's generator as it stands. A `posterior_function` is
called as it is and decides itself what to do with the seed; see
[`bayests_posterior`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md).

A sampler that cannot run raises its error rather than returning
something. The message names what about the input it could not work
with. Applied to a list of models – a 'modellist', an 'expandingwindow'
or, in bgvars, a 'gvecmodel' – that ends the run on the first
specification that fails, rather than leaving that one without a
posterior and carrying it into whatever reads the results.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for
VAR model restrictions. *Journal of Econometrics, 142*(1), 553–580.
[doi:10.1016/j.jeconom.2007.08.017](https://doi.org/10.1016/j.jeconom.2007.08.017)

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224–242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

## See also

[`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
describes the object this returns, element by element.

Other posterior simulation:
[`add_forecast_input.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md),
[`add_forecast_input.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md),
[`add_posterior_coefficients.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md),
[`add_posterior_forecasts.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md),
[`add_posterior_forecasts.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md),
[`add_posterior_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md),
[`add_posterior_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md),
[`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
[`bayests_files()`](https://franzmohr.github.io/bvartools/reference/bayests_files.md),
[`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md),
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`chain_diagnostics()`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)

## Examples

``` r

# Load data 
data("e6")
e6 <- e6 * 100

# Generate model
model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)
```
