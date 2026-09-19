# Thinning Posterior Draws

Thins the MCMC posterior draws in an object of class 'bvarmodel'.

## Usage

``` r
# S3 method for class 'bvarmodel'
thin(x, thin = 10, ...)
```

## Arguments

- x:

  an object of class 'bvarmodel'.

- thin:

  an integer specifying the thinning interval between successive values
  of posterior draws.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'bvarmodel'.

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

# Thinning
model <- thin(model, 2)
```
