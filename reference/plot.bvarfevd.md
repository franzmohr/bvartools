# Plotting Forecast Error Variance Decompositions of Bayesian Vector Autoregression

A plot function for objects of class "bvarfevd".

## Usage

``` r
# S3 method for class 'bvarfevd'
plot(x, max_groups = NULL, ...)
```

## Arguments

- x:

  an object of class "bvarfevd", usually, a result of a call to
  [`fevd`](https://franzmohr.github.io/bvartools/reference/fevd.md).

- max_groups:

  integer. Maximum number of variables shown in the plot. The
  `max_groups - 1` variables with the largest contributions across the
  whole horizon are kept and the contributions of the remaining
  variables are added up in a further bar segment named `"Other"`, so
  that the legend does not become too large. Default is `NULL`, so that
  a segment is shown for every variable of the decomposition.

- ...:

  further graphical parameters.

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100
e1 <- window(e1, end = c(1978, 4))

# Generate model data
model <- create_bvarmodel(e1, p = 2, deterministic = "const",
                          iterations = 100, burnin = 10)
# Chosen number of iterations and burnin should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws
object <- add_posterior_coefficients(model)

# Obtain FEVD
vd <- fevd(object, response = "cons")

# Plot
plot(vd)

```
