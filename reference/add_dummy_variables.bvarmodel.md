# Dummy Variables

Adds dummy variables to the deterministic terms of a model: impulse
dummies for single unusual periods, step dummies for level shifts, or
any other series that is known in advance.

## Usage

``` r
# S3 method for class 'bvarmodel'
add_dummy_variables(object, impulse = NULL, step = NULL, data = NULL, ...)

# S3 method for class 'bvecmodel'
add_dummy_variables(object, impulse = NULL, step = NULL, data = NULL, ...)

# S3 method for class 'expandingwindow'
add_dummy_variables(object, impulse = NULL, step = NULL, data = NULL, ...)

# S3 method for class 'modellist'
add_dummy_variables(object, ...)
```

## Arguments

- object:

  an object of class 'bvarmodel' or 'bvecmodel', usually the output of
  [`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  or
  [`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
  or a list of such models of class 'modellist' or 'expandingwindow'.

- impulse:

  a period, given the way [`ts`](https://rdrr.io/r/stats/ts.html) takes
  `start`, e.g. `c(2020, 2)`, or a list of such periods. Each gets a
  dummy variable that is one in that period and zero otherwise.

- step:

  a period or a list of periods as in `impulse`. Each gets a dummy
  variable that is zero before that period and one from it on.

- data:

  an optional time-series object of further dummy variables, or of any
  other series known in advance, with named columns at the frequency of
  the model. It has to cover the estimation sample. See 'Details'.

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with the dummy variables added to `data$train$x`
and its SUR form `data$train$z`, to `data$original$deterministic` and to
the names in `model$deterministic`, `model$n` increased by their number,
and `model$dummy_variables`, a data frame with one row per dummy
variable giving its `name`, its `type` (`"impulse"`, `"step"` or
`"data"`) and the `time` of its period, which is what a forecast
continues it by.

## Details

A dummy variable takes up what a model otherwise has no explanation for:
an impulse dummy the one period of a strike, a tax change or a pandemic,
a step dummy a lasting shift in the level of the series, such as a
change in how they are measured. Without one, a single extreme
observation pulls the coefficients of every equation towards itself and
inflates the estimated error variance for the whole sample.

The dummies are deterministic terms like the constant. They are appended
to the end of the deterministic terms of the regressors, named after
their type and period – `impulse.2020Q2`, `step.2008M09` – or after the
columns of `data`, and their coefficients take the prior of the
deterministic terms, `v_i_det` in
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).
**The coefficient of an impulse dummy is estimated from a single
observation**, so its prior matters more than that of the constant: a
tight one leaves most of the observation to the other coefficients,
which is what the dummy was meant to prevent.

**The function has to be called before
[`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md)**,
since it changes the number of coefficients the priors are given for. A
dummy variable that is zero in every period of the estimation sample, or
that the other regressors already span – a step dummy that is one in
every period is the constant again – is refused.

A forecast continues the dummies by their definition: an impulse dummy
is zero after its period and a step dummy stays one. The series in
`data` are continued with the values they were given, so for a forecast
they have to reach the last forecast period, or all deterministic terms
of the forecast have to be given in argument `deterministic` of
[`add_forecast_input`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md).

Of a list of models of class 'expandingwindow', each window gets the
dummies whose periods it reaches and leaves out the rest: a forecast
made at the end of a window that ends before a period could not have
known about what happened in it.
[`use_expanding_window`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
applied to a model with dummy variables does the same.

In a VEC model the dummies enter the equations of the differences, not
the cointegration term. **An impulse dummy there shifts the level of the
series for good**, since a one-off change in a difference is a lasting
change in the level. A one-off blip in the level is a dummy that is one
in its period and minus one in the next, which can be given in `data`.
[`vec_to_var`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
carries the dummies into the VAR representation.

## See also

Other model set-up:
[`add_initial_values.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md),
[`add_initial_values.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md),
[`add_prior_options()`](https://franzmohr.github.io/bvartools/reference/add_prior_options.md),
[`add_priors.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md),
[`add_priors.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md),
[`aggregation_weights()`](https://franzmohr.github.io/bvartools/reference/aggregation_weights.md),
[`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md),
[`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md),
[`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md),
[`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md),
[`use_expanding_window.bvarmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md),
[`use_expanding_window.bvecmodel()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# An impulse dummy for the first quarter of 1975 and a step dummy from 1979
model <- add_dummy_variables(model, impulse = c(1975, 1), step = c(1979, 1))
model[["model"]][["deterministic"]]
#> [1] "const"          "impulse.1975Q1" "step.1979Q1"   

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))
```
