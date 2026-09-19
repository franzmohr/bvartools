# Add Log-Likelihood

Forwards its elements to methods used to calculate posterior
log-likelihoods.

## Usage

``` r
# S3 method for class 'modellist'
add_posterior_loglik(object, ..., cores = 1)
```

## Arguments

- object:

  a list of class 'modellist', usually, the result of a call to
  [`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md).

- ...:

  further arguments passed to or from other methods.

- cores:

  the number of worker processes the models are simulated on. Defaults
  to 1, which simulates them one after the other in this session. See
  section 'Parallel simulation'.

## Value

The object in `object` with log-likelihood draws added to each of its
models, as described in
[`add_posterior_loglik.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md).

## Parallel simulation

With `cores` above 1 the models are simulated on a cluster of that many
worker processes, started with
[`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html) and stopped
when the call returns, but on no more workers than there are models.
Nested lists, such as the 'modellist' of 'expandingwindow' lists that
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
returns for several specifications, are shared out model by model.

Each worker runs with a single BLAS thread: the function sets
`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and
`VECLIB_MAXIMUM_THREADS` to 1 while the workers start, whatever this
session has set, and puts them back afterwards. An optimised BLAS would
otherwise start one thread per core in every worker, and the samplers do
not become faster with more threads. The thread count of this session is
not changed.

The workers use the library paths of this session and load the installed
version of the package, together with every package that provides a
method the models of the list need – dfmtools for its dynamic factor
models, for example. A function passed as `posterior_function`, such as
the one
[`bayests_posterior`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md)
returns, is sent to the workers along with the variables of the
environment it was created in, and must not rely on anything else of
this session.

Coefficient draws do not depend on the number of workers. Every model is
simulated with its own seed, see
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md),
and a model without one is given one from R's random number generator
before it is sent off. The draws equal those of a simulation on one core
in a session whose BLAS runs on one thread, as R's reference BLAS always
does. An optimised BLAS such as OpenBLAS rounds some results differently
in the last digits on one thread than on several, and a Markov chain
carries such a difference forward, so a session running it on several
threads gets different draws from the same posterior. Draws that are not
seeded per model, such as those of
[`add_posterior_forecasts`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md),
come from independent random number streams of the workers, set up with
[`clusterSetRNGStream`](https://rdrr.io/r/parallel/RngStream.html) from
R's generator. They are reproduced by
[`set.seed`](https://rdrr.io/r/base/Random.html) for a given number of
workers, but differ from those of a simulation on one core.

Starting the workers takes a moment, so a cluster pays off for lists of
models that take longer than that to simulate.

## Examples

``` r

# Load data 
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model
model <- create_bvarmodel(e1, p = 1:2, deterministic = 2,
                          iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
object <- add_posterior_coefficients(model)

# Add log-likelihoods
object <- add_posterior_loglik(object)
```
