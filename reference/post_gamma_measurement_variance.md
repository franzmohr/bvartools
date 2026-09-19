# Posterior Draws of Error Variances

Produces a draw of the constant diagonal error variance matrix of the
measurement equation of a state space model using an inverse gamma
posterior density.

## Usage

``` r
post_gamma_measurement_variance(u, shape_prior, rate_prior, inverse)
```

## Arguments

- u:

  a \\KT \times 1\\ vector of errors.

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

For a model with measurement equation \$\$y_t = Z\_{t} a_t + u_t\$\$
with \\u_t \sim N(0, \Sigma\_{u})\\ the function produces a draw of the
constant diagonal error variance matrix \\\Sigma_u\\.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

## Examples

``` r

k <- 10 # Number of endogenous variables
tt <- 1000 # Number of observations

set.seed(1234) # Set RNG seed

# Generate artificial error series with N(0, 1)
u <- matrix(rnorm(k * tt))

# Define priors
shape_prior <- matrix(1, k)
rate_prior <- matrix(.0001, k)

# Obtain posterior draw
post_gamma_measurement_variance(u, shape_prior, rate_prior, inverse = FALSE)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                                                        
#>  [1,] 1.080032 .       .         .         .         .         .       
#>  [2,] .        1.03679 .         .         .         .         .       
#>  [3,] .        .       0.9076317 .         .         .         .       
#>  [4,] .        .       .         0.9669023 .         .         .       
#>  [5,] .        .       .         .         0.8980991 .         .       
#>  [6,] .        .       .         .         .         0.9443887 .       
#>  [7,] .        .       .         .         .         .         0.950957
#>  [8,] .        .       .         .         .         .         .       
#>  [9,] .        .       .         .         .         .         .       
#> [10,] .        .       .         .         .         .         .       
#>                                   
#>  [1,] .         .         .       
#>  [2,] .         .         .       
#>  [3,] .         .         .       
#>  [4,] .         .         .       
#>  [5,] .         .         .       
#>  [6,] .         .         .       
#>  [7,] .         .         .       
#>  [8,] 0.9798014 .         .       
#>  [9,] .         0.9263457 .       
#> [10,] .         .         1.042936
```
