# Plotting Draws of a Bayesian VAR Model

A plot function for objects of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
plot(x, ci = 0.95, type = "hist", show_zero_y = TRUE, max_cols = 6, ...)
```

## Arguments

- x:

  an object of class 'bvarmodel'.

- ci:

  interval used to calculate credible bands for time-varying parameters.

- type:

  either `"hist"` (default) for histograms, `"trace"` for a trace plot
  or `"boxplot"` for a boxplot. Only used for parameter draws of
  constant coefficients.

- show_zero_y:

  if `TRUE` (default), a horizontal line with y = 0 is added to the
  plot. Only used for time varying parameters.

- max_cols:

  an integer of the maximum number of regressors per figure. A block
  with more regressors than this is drawn as several figures of nearly
  equal width. Defaults to 6.

- ...:

  further graphical parameters.

## Value

A plot per block of coefficients.

## Details

The function draws one figure per block of coefficients – the lags of
the endogenous variables, the exogenous variables, the deterministic
terms, the contemporaneous endogenous variables of a structural model,
and the covariance matrix of the error term – instead of one figure for
the whole model. A model with many regressors would otherwise produce
panels too small to read.

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
model <- add_posterior_coefficients(model)

# Plot
plot(model, type = "hist")



plot(model, type = "trace")



plot(model, type = "boxplot")




```
