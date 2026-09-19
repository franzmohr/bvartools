# SUR Matrix Transformation

Transforms a dense matrix of dimensions \\KT \times M\\ into a sparse
block diagonal matrix of dimensions \\KT \times MT\\.

## Usage

``` r
sur_const_to_tvp(z, k, tt)
```

## Arguments

- z:

  a \\KT \times M\\ matrix.

- k:

  integer of the number of endogenous variables.

- tt:

  integer of the number of observations.

## Value

A sparse block diagonal matrix of dimensions \\KT \times MT\\.

## Examples

``` r

# Specify the dimensions of the dense matrix
k <- 2
tt <- 5
m <- 3

# Generate artificial data
z <- matrix(NA, k * tt, m)
for (i in 1:tt) {
  z[(i - 1) * k + 1:k, ] <- i
}

# Perform transformation
sur_const_to_tvp(z, k, tt)
#> 10 x 15 sparse Matrix of class "dgCMatrix"
#>                                    
#>  [1,] 1 1 1 . . . . . . . . . . . .
#>  [2,] 1 1 1 . . . . . . . . . . . .
#>  [3,] . . . 2 2 2 . . . . . . . . .
#>  [4,] . . . 2 2 2 . . . . . . . . .
#>  [5,] . . . . . . 3 3 3 . . . . . .
#>  [6,] . . . . . . 3 3 3 . . . . . .
#>  [7,] . . . . . . . . . 4 4 4 . . .
#>  [8,] . . . . . . . . . 4 4 4 . . .
#>  [9,] . . . . . . . . . . . . 5 5 5
#> [10,] . . . . . . . . . . . . 5 5 5

```
