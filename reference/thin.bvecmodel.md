# Thinning Posterior Draws

Thins the MCMC posterior draws in an object of class 'bvecmodel'.

## Usage

``` r
# S3 method for class 'bvecmodel'
thin(x, thin = 10, ...)
```

## Arguments

- x:

  an object of class 'bvecmodel'.

- thin:

  an integer specifying the thinning interval between successive values
  of posterior draws.

- ...:

  further arguments passed to or from other methods.

## Value

An object of class 'bvecmodel' with every block of draws thinned.

## Examples

``` r
# Load data
data("e6")

# Create model
model <- create_bvecmodel(e6, p = 2, r = 1,
                         const = "unrestricted",
                          iterations = 20, burnin = 10)
# Number of iterations and burnin should be much higher.

# Add priors
model <- add_priors(model,
                    coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))

# Add initial values
model <- add_initial_values(model)

# Obtain posterior draws 
model <- add_posterior_coefficients(model)

# Thinning
model <- thin(model, 2)
```
