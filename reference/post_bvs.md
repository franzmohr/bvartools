# Bayesian Variable Selection

`post_bvs` employs Bayesian variable selection as proposed by Korobilis
(2013) to produce a vector of inclusion parameters for the coefficient
matrix of a VAR model.

## Usage

``` r
post_bvs(y, z, a, k, m, lambda, sigma_i, prob_prior, include = NULL)
```

## Arguments

- y:

  a \\KT \times 1\\ vector of the endogenous variables.

- z:

  a \\KT \times M\\ matrix of explanatory variables.

- a:

  an M-dimensional vector of parameter draws. If time varying parameters
  are used, an \\M \times T\\ coefficient matrix can be provided.

- k:

  integer of the number of endogenous variables.

- m:

  integer of the number of M

- lambda:

  an \\M \times M\\ inclusion matrix that should be updated.

- sigma_i:

  a sparse \\KT \times KT\\ block diagonal matrix containing the inverse
  variance-covariance matrix.

- prob_prior:

  an M-dimensional vector of prior inclusion probabilities.

- include:

  an integer vector specifying the positions of variables, which should
  be included in the BVS algorithm. If `NULL` (default), BVS will be
  applied to all variables.

## Value

A matrix of inclusion parameters on its diagonal.

## Details

The function employs Bayesian variable selection as proposed by
Korobilis (2013) to produce a vector of inclusion parameters, which are
the diagonal elements of the inclusion matrix \\\Lambda\\ for the VAR
model \$\$y_t = Z_t \Lambda a_t + u_t,\$\$ where \\u_t \sim N(0,
\Sigma\_{t})\\. \\y_t\\ is a K-dimensional vector of endogenous
variables and \\Z_t = x_t^{\prime} \otimes I_K\\ is a \\K \times M\\
matrix of regressors with \\x_t\\ as a vector of regressors.

## References

Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
*Journal of Applied Econometrics, 28*(2), 204–230.
[doi:10.1002/jae.1271](https://doi.org/10.1002/jae.1271)

## Examples

``` r

# Load data
data("e1")
e1 <- diff(log(e1)) * 100

# Generate model input, which uses BVS as variable selection algorithm
object <- create_bvarmodel(data = e1, p = 2, deterministic = "const",
                           varsel = "bvs")

# Add prior specifications, including the prior inclusion probabilities
object <- add_priors(object,
                     coef = list(v_i = 1, v_i_det = 1 / 10),
                     sigma = list(df = "k", scale = 1),
                     varsel = list(inprior = .1))

# Add initial values
object <- add_initial_values(object)

# Obtain data, initial values and priors
y <- matrix(t(object[["data"]][["train"]][["y"]]))
z <- object[["data"]][["train"]][["z"]] # Argument 'z' is taken dense
k <- object[["model"]][["k"]]
tt <- nrow(object[["data"]][["train"]][["y"]])
m <- ncol(z)
a <- object[["initial"]][["a"]]
prob_prior <- object[["priors"]][["a"]][["inprior"]]

# Arguments 'lambda' and 'sigma_i' have to be sparse
lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)

# Initial value of the inverse error covariance matrix
u <- matrix(y - z %*% a, k)
sigma_i <- Matrix::Matrix(kronecker(diag(1, tt), solve(tcrossprod(u) / tt)),
                          sparse = TRUE)

# Draw inclusion parameters
post_bvs(y, z, a, k, m, lambda, sigma_i, prob_prior)
#> 21 x 21 sparse Matrix of class "dgCMatrix"
#>                                                
#>  [1,] 1 . . . . . . . . . . . . . . . . . . . .
#>  [2,] . 1 . . . . . . . . . . . . . . . . . . .
#>  [3,] . . . . . . . . . . . . . . . . . . . . .
#>  [4,] . . . 1 . . . . . . . . . . . . . . . . .
#>  [5,] . . . . 1 . . . . . . . . . . . . . . . .
#>  [6,] . . . . . 1 . . . . . . . . . . . . . . .
#>  [7,] . . . . . . 1 . . . . . . . . . . . . . .
#>  [8,] . . . . . . . 1 . . . . . . . . . . . . .
#>  [9,] . . . . . . . . 1 . . . . . . . . . . . .
#> [10,] . . . . . . . . . . . . . . . . . . . . .
#> [11,] . . . . . . . . . . 1 . . . . . . . . . .
#> [12,] . . . . . . . . . . . 1 . . . . . . . . .
#> [13,] . . . . . . . . . . . . . . . . . . . . .
#> [14,] . . . . . . . . . . . . . . . . . . . . .
#> [15,] . . . . . . . . . . . . . . 1 . . . . . .
#> [16,] . . . . . . . . . . . . . . . 1 . . . . .
#> [17,] . . . . . . . . . . . . . . . . 1 . . . .
#> [18,] . . . . . . . . . . . . . . . . . 1 . . .
#> [19,] . . . . . . . . . . . . . . . . . . 1 . .
#> [20,] . . . . . . . . . . . . . . . . . . . 1 .
#> [21,] . . . . . . . . . . . . . . . . . . . . 1
```
