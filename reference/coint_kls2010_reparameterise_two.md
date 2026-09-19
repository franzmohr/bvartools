# Cointegration Reparameterisation

Performs the second transformation of the loading and cointegration
matrix as proposed in Koop et al. (2010).

## Usage

``` r
coint_kls2010_reparameterise_two(alpha, beta)
```

## Arguments

- alpha:

  a \\K \times r\\ matrix.

- beta:

  an \\M \times r\\ matrix.

## Value

A list of two matrices:

- alpha:

  Loading matrix \\A\\.

- beta:

  Cointegration matrix \\B\\ (semiorthogonal).

## Details

The function performs two transformations:

- \\A = \alpha (\beta^{\prime} \beta)^{1/2}\\

- \\B = \beta (\beta^{\prime} \beta)^{-1/2}\\

## References

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224-242.
[doi:10.1080/07474930903382208](https://doi.org/10.1080/07474930903382208)

## Examples

``` r

# Generate input data
alpha <- matrix(c(-0.07, 0.17), 2)
beta <- matrix(c(1, -4), 2)

# Reparameterise
coint_kls2010_reparameterise_two(alpha, beta)
#> $alpha
#>            [,1]
#> [1,] -0.2886174
#> [2,]  0.7009280
#> 
#> $beta
#>            [,1]
#> [1,]  0.2425356
#> [2,] -0.9701425
#> 
```
