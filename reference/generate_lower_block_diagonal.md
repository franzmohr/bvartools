# Posterior Data Preparation

Generates a lower triangular block matrix with ones on the main
diagonal, where the off-diagonal elements are blocks.

## Usage

``` r
generate_lower_block_diagonal(a, k, tt)
```

## Arguments

- a:

  \\K^2p\\-dimensional vector of coefficients. See 'Details'.

- k:

  integer of the number of columns per block.

- tt:

  integer specifying

## Value

A sparse matrix.

## Details

For the \\K \times Kp\\ matrix A, where \\a = vec(A)\\, with \\A =
\left\[A_1, A_2, ..., A_p\right\]\\ the function constructs the
following sparse \\KT \times KT\\ diagonal block matrix:
\$\$\begin{bmatrix} I\_{K} & 0 & 0 & 0 & \dots & 0 \\-A\_{1} & I\_{K} &
0 & 0 & \dots & 0 \\-A\_{2} & -A\_{1}& I\_{K} & 0 & \dots & 0 \\ 0 &
-A\_{2}& -A\_{1}& I\_{K} & \dots & 0 \\ \vdots & \ddots& \ddots& \ddots&
\ddots& 0 \\ 0 & \dots & 0 & -A\_{2}& -A\_{1}& I\_{K} \end{bmatrix}.\$\$

## Examples

``` r
a <- matrix(1:8)
generate_lower_block_diagonal(a, 2, 5)
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                  
#>  [1,]  1  .  .  .  .  .  .  . . .
#>  [2,]  .  1  .  .  .  .  .  . . .
#>  [3,] -1 -3  1  .  .  .  .  . . .
#>  [4,] -2 -4  .  1  .  .  .  . . .
#>  [5,] -5 -7 -1 -3  1  .  .  . . .
#>  [6,] -6 -8 -2 -4  .  1  .  . . .
#>  [7,]  .  . -5 -7 -1 -3  1  . . .
#>  [8,]  .  . -6 -8 -2 -4  .  1 . .
#>  [9,]  .  .  .  . -5 -7 -1 -3 1 .
#> [10,]  .  .  .  . -6 -8 -2 -4 . 1
```
