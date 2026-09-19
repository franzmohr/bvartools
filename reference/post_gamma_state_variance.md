# Posterior Draws of Error Variances

Produces a draw of the constant diagonal error variance matrix of the
state equation of a state space model using an inverse gamma posterior
density.

## Usage

``` r
post_gamma_state_variance(a, a_init, shape_prior, rate_prior, inverse)
```

## Arguments

- a:

  a \\KT \times 1\\ vector of time varying parameter draws.

- a_init:

  a \\K \times 1\\ vector of initial states.

- shape_prior:

  a \\K \times 1\\ vector of prior shape parameters.

- rate_prior:

  a \\K \times 1\\ vector of prior rate parameters.

- inverse:

  logical. If `TRUE`, the function returns the precision matrix, i.e.
  the inverse of the variance matrix. Defaults to `FALSE`.

## Value

A matrix.

## Details

For the state space model with state equation \$\$a_t = a\_{t-1} +
v_t\$\$ and measurement equation \$\$y_t = Z\_{t} a\_{t} + u\_{t}\$\$
with \\v_t \sim N(0, \Sigma^{v})\\ and \\u_t \sim N(0, \Sigma^{u})\\ the
function produces a draw of the constant diagonal error variances matrix
of the state equation \\\Sigma^{v}\\.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

## Examples

``` r
k <- 10 # Number of artificial coefficients
tt <- 1000 # Number of observations

set.seed(1234) # Set RNG seed

# Generate artificial data according to a random walk
a <- matrix(rnorm(k), k, tt + 1)
for (i in 2:(tt + 1)) {
  a[, i] <- a[, i - 1] + rnorm(k, 0, sqrt(1 / 100))
}

a_init <- matrix(a[, 1]) # Define initial state
a <- matrix(a[, -1]) # Drop initial state from main sample and make vector

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
post_gamma_state_variance(a, a_init, shape_prior, rate_prior, inverse = FALSE)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                          
#>  [1,] 0.0103381 .          .           .           .           .         
#>  [2,] .         0.01019774 .           .           .           .         
#>  [3,] .         .          0.009851695 .           .           .         
#>  [4,] .         .          .           0.009406544 .           .         
#>  [5,] .         .          .           .           0.009784046 .         
#>  [6,] .         .          .           .           .           0.01023659
#>  [7,] .         .          .           .           .           .         
#>  [8,] .         .          .           .           .           .         
#>  [9,] .         .          .           .           .           .         
#> [10,] .         .          .           .           .           .         
#>                                                     
#>  [1,] .           .           .           .         
#>  [2,] .           .           .           .         
#>  [3,] .           .           .           .         
#>  [4,] .           .           .           .         
#>  [5,] .           .           .           .         
#>  [6,] .           .           .           .         
#>  [7,] 0.009754256 .           .           .         
#>  [8,] .           0.009246524 .           .         
#>  [9,] .           .           0.008473591 .         
#> [10,] .           .           .           0.00947631
```
