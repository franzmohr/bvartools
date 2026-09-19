# Covariance: Data Preparation

Convenience function, which generates the input data for posterior
simulation of covariance parameters.

## Usage

``` r
covar_prepare_data(y, omega_i, k, tt, tvp)
```

## Arguments

- y:

  a \\KT \times 1\\ vector of input data.

- omega_i:

  a \\K \times K\\ or \\KT \times KT\\ matrix of error variances. The
  matrix must be sparse.

- k:

  an integer of the number of endogenous variables.

- tt:

  an integer of the number of observations.

- tvp:

  logical indicating if the SUR matrix with the values of regressors
  should be prepared for the estimation of constant or time varying
  parameters.

## Value

A list with three elements:

- y:

  The prepared vector of endogenous variables.

- z:

  The prepared matrix of regressors.

- omega_i:

  The prepared diagonal matrix of measurement error variances.

All matrices are returned as sparse matrices.

## Details

For the model \$\$y_t = Z\_{t} a_t + u_t\$\$ with \\u_t \sim N(0, \Psi
\Omega\_{t} \Psi^{\prime})\\ and \\\Omega\_{t}\\ as a diagonal matrix of
error variances, the function produces the input data for the posterior
simulation of the lower triangular covariance coefficients of \\\Psi\\
as presented in Primiceri (2005).

## References

Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). *Bayesian
econometric methods* (2nd ed.). Cambridge: Cambridge University Press.

Primiceri, G. E. (2005). Time varying structural vector autoregressions
and monetary policy. *The Review of Economic Studies 72*(3), 821–852.
[doi:10.1111/j.1467-937X.2005.00353.x](https://doi.org/10.1111/j.1467-937X.2005.00353.x)

## Examples

``` r

k <- 3
tt <- 4
u <- matrix(1:(k * tt))
  
# Generate simple variance matrix
omega_i <- Matrix::Matrix(diag(1:k, k))
# Generate block diagonal variance matrix
tv_omega_i <- Matrix::Matrix(0, k * tt, k * tt)
for (i in 1:tt) {
  tv_omega_i[(i - 1) * k + 1:k, (i - 1) * k + 1:k] <- omega_i
}

# Constant error variances

# Constant coefficients
covar_prepare_data(u, omega_i, k, tt, FALSE)
#> $y
#>      [,1]
#> [1,]    2
#> [2,]    3
#> [3,]    5
#> [4,]    6
#> [5,]    8
#> [6,]    9
#> [7,]   11
#> [8,]   12
#> 
#> $z
#> 8 x 3 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,]  -1   .   .
#> [2,]   .  -1  -2
#> [3,]  -4   .   .
#> [4,]   .  -4  -5
#> [5,]  -7   .   .
#> [6,]   .  -7  -8
#> [7,] -10   .   .
#> [8,]   . -10 -11
#> 
#> $omega_i
#> 8 x 8 sparse Matrix of class "dgCMatrix"
#>                     
#> [1,] 2 . . . . . . .
#> [2,] . 3 . . . . . .
#> [3,] . . 2 . . . . .
#> [4,] . . . 3 . . . .
#> [5,] . . . . 2 . . .
#> [6,] . . . . . 3 . .
#> [7,] . . . . . . 2 .
#> [8,] . . . . . . . 3
#> 
# Time varying coefficients
covar_prepare_data(u, omega_i, k, tt, TRUE)
#> $y
#>      [,1]
#> [1,]    2
#> [2,]    3
#> [3,]    5
#> [4,]    6
#> [5,]    8
#> [6,]    9
#> [7,]   11
#> [8,]   12
#> 
#> $z
#> 8 x 12 sparse Matrix of class "dgCMatrix"
#>                                            
#> [1,] -1  .  .  .  .  .  .  .  .   .   .   .
#> [2,]  . -1 -2  .  .  .  .  .  .   .   .   .
#> [3,]  .  .  . -4  .  .  .  .  .   .   .   .
#> [4,]  .  .  .  . -4 -5  .  .  .   .   .   .
#> [5,]  .  .  .  .  .  . -7  .  .   .   .   .
#> [6,]  .  .  .  .  .  .  . -7 -8   .   .   .
#> [7,]  .  .  .  .  .  .  .  .  . -10   .   .
#> [8,]  .  .  .  .  .  .  .  .  .   . -10 -11
#> 
#> $omega_i
#> 8 x 8 sparse Matrix of class "dgCMatrix"
#>                     
#> [1,] 2 . . . . . . .
#> [2,] . 3 . . . . . .
#> [3,] . . 2 . . . . .
#> [4,] . . . 3 . . . .
#> [5,] . . . . 2 . . .
#> [6,] . . . . . 3 . .
#> [7,] . . . . . . 2 .
#> [8,] . . . . . . . 3
#> 


# Time varying error variances

# Constant coefficients
covar_prepare_data(u, tv_omega_i, k, tt, FALSE)
#> $y
#>      [,1]
#> [1,]    2
#> [2,]    3
#> [3,]    5
#> [4,]    6
#> [5,]    8
#> [6,]    9
#> [7,]   11
#> [8,]   12
#> 
#> $z
#> 8 x 3 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,]  -1   .   .
#> [2,]   .  -1  -2
#> [3,]  -4   .   .
#> [4,]   .  -4  -5
#> [5,]  -7   .   .
#> [6,]   .  -7  -8
#> [7,] -10   .   .
#> [8,]   . -10 -11
#> 
#> $omega_i
#> 8 x 8 sparse Matrix of class "dgCMatrix"
#>                     
#> [1,] 2 . . . . . . .
#> [2,] . 3 . . . . . .
#> [3,] . . 2 . . . . .
#> [4,] . . . 3 . . . .
#> [5,] . . . . 2 . . .
#> [6,] . . . . . 3 . .
#> [7,] . . . . . . 2 .
#> [8,] . . . . . . . 3
#> 
# Time varying coefficients
covar_prepare_data(u, tv_omega_i, k, tt, TRUE)
#> $y
#>      [,1]
#> [1,]    2
#> [2,]    3
#> [3,]    5
#> [4,]    6
#> [5,]    8
#> [6,]    9
#> [7,]   11
#> [8,]   12
#> 
#> $z
#> 8 x 12 sparse Matrix of class "dgCMatrix"
#>                                            
#> [1,] -1  .  .  .  .  .  .  .  .   .   .   .
#> [2,]  . -1 -2  .  .  .  .  .  .   .   .   .
#> [3,]  .  .  . -4  .  .  .  .  .   .   .   .
#> [4,]  .  .  .  . -4 -5  .  .  .   .   .   .
#> [5,]  .  .  .  .  .  . -7  .  .   .   .   .
#> [6,]  .  .  .  .  .  .  . -7 -8   .   .   .
#> [7,]  .  .  .  .  .  .  .  .  . -10   .   .
#> [8,]  .  .  .  .  .  .  .  .  .   . -10 -11
#> 
#> $omega_i
#> 8 x 8 sparse Matrix of class "dgCMatrix"
#>                     
#> [1,] 2 . . . . . . .
#> [2,] . 3 . . . . . .
#> [3,] . . 2 . . . . .
#> [4,] . . . 3 . . . .
#> [5,] . . . . 2 . . .
#> [6,] . . . . . 3 . .
#> [7,] . . . . . . 2 .
#> [8,] . . . . . . . 3
#> 
```
