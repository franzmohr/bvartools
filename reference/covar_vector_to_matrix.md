# Covariance: Vector to Matrix

Convenience function, which takes the vector of draws of lower
triangular covariance coefficients and transforms it into a matrix with
ones on the main diagonal. In case of time varying parameters the
resulting matrix will be block diagonal.

## Usage

``` r
covar_vector_to_matrix(psi, k, tt)
```

## Arguments

- psi:

  a \\K (K - 1) / 2 \times 1\\ or \\T K (K - 1) / 2 \times 1\\ vector of
  input data.

- k:

  the number \\K\\ of endogenous variables.

- tt:

  the number \\T\\ of observations.

## Value

A sparse, block diagonal matrix.

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

Primiceri, G. E. (2005). Time varying structural vector autoregressions
and monetary policy. *The Review of Economic Studies 72*(3), 821–852.
[doi:10.1111/j.1467-937X.2005.00353.x](https://doi.org/10.1111/j.1467-937X.2005.00353.x)

## Examples

``` r

# Create artificial data
k <- 5
tt <- 4
n_covar <- (k - 1) * k / 2

# Constant parameters
psi <- matrix(1:(n_covar))
covar_vector_to_matrix(psi, k, tt)
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                
#> [1,] 1 . .  . .
#> [2,] 1 1 .  . .
#> [3,] 2 3 1  . .
#> [4,] 4 5 6  1 .
#> [5,] 7 8 9 10 1

# Time varying parameters
psi <- matrix(1:(n_covar * tt))
covar_vector_to_matrix(psi, k, tt)
#> 20 x 20 sparse Matrix of class "dgCMatrix"
#>                                                           
#>  [1,] 1 . .  . .  .  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [2,] 1 1 .  . .  .  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [3,] 2 3 1  . .  .  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [4,] 4 5 6  1 .  .  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [5,] 7 8 9 10 1  .  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [6,] . . .  . .  1  .  .  . .  .  .  .  . .  .  .  .  . .
#>  [7,] . . .  . . 11  1  .  . .  .  .  .  . .  .  .  .  . .
#>  [8,] . . .  . . 12 13  1  . .  .  .  .  . .  .  .  .  . .
#>  [9,] . . .  . . 14 15 16  1 .  .  .  .  . .  .  .  .  . .
#> [10,] . . .  . . 17 18 19 20 1  .  .  .  . .  .  .  .  . .
#> [11,] . . .  . .  .  .  .  . .  1  .  .  . .  .  .  .  . .
#> [12,] . . .  . .  .  .  .  . . 21  1  .  . .  .  .  .  . .
#> [13,] . . .  . .  .  .  .  . . 22 23  1  . .  .  .  .  . .
#> [14,] . . .  . .  .  .  .  . . 24 25 26  1 .  .  .  .  . .
#> [15,] . . .  . .  .  .  .  . . 27 28 29 30 1  .  .  .  . .
#> [16,] . . .  . .  .  .  .  . .  .  .  .  . .  1  .  .  . .
#> [17,] . . .  . .  .  .  .  . .  .  .  .  . . 31  1  .  . .
#> [18,] . . .  . .  .  .  .  . .  .  .  .  . . 32 33  1  . .
#> [19,] . . .  . .  .  .  .  . .  .  .  .  . . 34 35 36  1 .
#> [20,] . . .  . .  .  .  .  . .  .  .  .  . . 37 38 39 40 1
```
