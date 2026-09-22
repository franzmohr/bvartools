# Posterior Simulation with the BayesTS Executable

Returns a function that simulates the posterior of a model with the
standalone BayesTS executable rather than with the samplers compiled
into this package.

## Usage

``` r
bayests_posterior(executable = NULL, library_path = NULL, scratch = tempdir())
```

## Arguments

- executable:

  the path to the BayesTS executable. If `NULL` (default), option
  `bvartools.bayests_executable` or, if that is not set, environment
  variable `BAYESTS_EXECUTABLE`.

- library_path:

  a character vector of directories the executable needs on its search
  path for its runtime libraries, put in front of `PATH` while it runs.
  A packaged BayesTS carries its libraries and does not need it; an
  executable from a build tree may.

- scratch:

  the directory where the model files are written while BayesTS draws
  into them. They are removed afterwards.

## Value

A function of one argument, a model of class 'bvarmodel' or 'bvecmodel',
that returns that model with its posterior draws.

## Details

The result is meant for argument `posterior_function` of
[`add_posterior_coefficients`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md):

    model <- add_posterior_coefficients(model,
      posterior_function = bayests_posterior())

It writes the model with
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md),
runs `bayests posterior` on the file without its log-likelihood and
forecast steps, reads the draws back with
[`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
and adds them to the model as element `posterior`. Nothing else about
the model is changed. The draws have the elements, dimensions and
thinning of those of the internal samplers, whose C++ code BayesTS
shares; a BayesTS build linked against an optimised BLAS draws faster.

BayesTS draws with the seed in `object$model$seed`, which
[`add_initial_values`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
sets and
[`add_seed`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
replaces. A model without one is given a seed drawn from R's random
number generator, since BayesTS would otherwise start every model from
the same state of its own generator. The seed used is kept in the
returned model.

The executable runs single threaded unless `OMP_NUM_THREADS` says
otherwise, which suits simulating one model per worker of a cluster.

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
[`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md),
[`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md),
[`chain_diagnostics()`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md),
[`predict.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)
