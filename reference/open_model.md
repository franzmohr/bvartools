# Open a Model Stored in an HDF5 File

Opens a model written with
[`write_to_hdf5`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
and returns a handle to it, which the analysis functions work through a
piece of the chain at a time rather than reading every draw into the
session.

## Usage

``` r
open_model(filename, group = "")

# S3 method for class 'bvarfile'
print(x, ...)
```

## Arguments

- filename:

  path to an HDF5 file holding a model.

- group:

  the group the model's tree hangs under inside its file, as
  [`read_model_from_hdf5`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
  takes it. Defaults to `""`, the root.

- x:

  an object of class 'bvarfile'.

- ...:

  further arguments passed to or from other methods.

## Value

A list of class 'bvarfile' with the elements `filename`, `group`,
`model`, `data`, `priors` and `draws`, the number of draws in the file.

## Details

The draws of a model are what it is large in. A model with time varying
coefficients keeps a path per draw, and a global model solved from many
sub-models keeps a square matrix per lag and draw, so a chain that a
posterior summary is comfortable with is a chain an analysis cannot hold
beside everything else it needs.

A handle carries what a model is – its specification, its data and its
priors – and the length of its chain, but none of the draws. The methods
for it read the draws in pieces and keep only what they are asked for:
the responses of an impulse response, the shares of a variance
decomposition. What they return is what the same call on the model in
memory returns.

## See also

Other model comparison:
[`add_forecast_errors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md),
[`add_forecast_errors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md),
[`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md),
[`add_predictive_loglik.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md),
[`add_predictive_loglik.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md),
[`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md),
[`align_model_obs.modellist()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md),
[`analysis_of_stored_models`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md),
[`choose_best_model.selcritlist()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md),
[`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md),
[`folder_steps`](https://franzmohr.github.io/bvartools/reference/folder_steps.md),
[`map_draws()`](https://franzmohr.github.io/bvartools/reference/map_draws.md),
[`map_models()`](https://franzmohr.github.io/bvartools/reference/map_models.md),
[`open_models()`](https://franzmohr.github.io/bvartools/reference/open_models.md),
[`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md),
[`selection_criteria.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.default()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md),
[`selection_criteria.modellist()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 10)
# Number of iterations and burn-in should be much higher.

model <- add_priors(model,
                    coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = .0001))
model <- add_posterior_coefficients(add_initial_values(model))

file <- file.path(tempdir(), "bvartools-example-model.h5")
unlink(file)
write_to_hdf5(model, filename = file)

stored <- open_model(file)
stored
#> Model in /tmp/RtmpcU1dnC/bvartools-example-model.h5 
#> 100 draws of a VarNormalWishart model of 3 variables

# The analysis reads the draws in pieces
ir <- irf(stored, impulse = "invest", response = "income", n_ahead = 5,
          chunk = 25)
```
