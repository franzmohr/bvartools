# Add Initial Values of an MCMC Chain

Adds initial values to a list of models by passing each element to the
respective method.

## Usage

``` r
# S3 method for class 'modellist'
add_initial_values(object, ...)
```

## Arguments

- object:

  a list of class 'modellist'
  [`add_priors`](https://franzmohr.github.io/bvartools/reference/add_priors.md).

- ...:

  further arguments passed to or from other methods.

## Value

The object in `object` with initial values added to each of its models,
as described in
[`add_initial_values.bvarmodel`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md).

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Create model
model <- create_bvarmodel(e1, p = 0:2, deterministic = "const",
                          iterations = 10, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add prior specifications
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)
```
