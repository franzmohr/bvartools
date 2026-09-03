test_that("covar_vector_to_matrix builds a unit lower triangular matrix", {
  skip_without_matrix()

  # Three covariance coefficients fill the strict lower triangle of a 3 x 3
  # matrix, column by column.
  result <- as.matrix(covar_vector_to_matrix(matrix(1:3), 3, 2))

  expect_identical(dim(result), c(3L, 3L))
  expect_equal(result, matrix(c(1, 1, 2,
                                0, 1, 3,
                                0, 0, 1), 3))
})

test_that("time varying covariance coefficients become block diagonal", {
  skip_without_matrix()

  result <- as.matrix(covar_vector_to_matrix(matrix(1:6), 3, 2))

  expect_identical(dim(result), c(6L, 6L))
  # One block per period, each built from its own three coefficients.
  expect_equal(result[1:3, 1:3],
               as.matrix(covar_vector_to_matrix(matrix(1:3), 3, 1)))
  expect_equal(result[4:6, 4:6],
               as.matrix(covar_vector_to_matrix(matrix(4:6), 3, 1)))
  expect_true(all(result[1:3, 4:6] == 0))
  expect_true(all(result[4:6, 1:3] == 0))
})

test_that("generate_lower_block_diagonal places the negative lag blocks", {
  skip_without_matrix()

  # a = vec([A_1, A_2]) for k = 2.
  result <- as.matrix(generate_lower_block_diagonal(matrix(1:8), 2, 4))
  a1 <- matrix(1:4, 2)
  a2 <- matrix(5:8, 2)

  expect_identical(dim(result), c(8L, 8L))
  expect_equal(result[1:2, 1:2], diag(2))
  expect_equal(result[3:4, 1:2], -a1)
  expect_equal(result[5:6, 1:2], -a2)
  expect_equal(result[5:6, 3:4], -a1)
  expect_equal(result[7:8, 3:4], -a2)
  # The first block row has no lags to subtract yet.
  expect_true(all(result[1:2, 3:8] == 0))
  # A_2 does not reach further back than two periods.
  expect_true(all(result[7:8, 1:2] == 0))
})

test_that("sur_const_to_tvp spreads the regressors over the diagonal", {
  skip_without_matrix()

  z <- matrix(1:8, 4, 2)
  result <- as.matrix(sur_const_to_tvp(z, 2, 2))

  expect_identical(dim(result), c(4L, 4L))
  expect_equal(result[1:2, 1:2], z[1:2, ])
  expect_equal(result[3:4, 3:4], z[3:4, ])
  expect_true(all(result[1:2, 3:4] == 0))
  expect_true(all(result[3:4, 1:2] == 0))
})

test_that("covar_prepare_data returns the pieces of the covariance regression", {
  skip_without_matrix()

  k <- 3
  tt <- 4
  u <- matrix(1:(k * tt))
  omega_i <- Matrix::Matrix(diag(1:k, k))

  prepared <- covar_prepare_data(u, omega_i, k, tt, FALSE)

  expect_named(prepared, c("y", "z", "omega_i"))
  # The first variable of each period carries no covariance regressor, so it
  # drops out of the regression.
  expect_identical(nrow(prepared[["y"]]), as.integer((k - 1) * tt))
  expect_identical(ncol(as.matrix(prepared[["z"]])),
                   as.integer(k * (k - 1) / 2))
  expect_identical(dim(as.matrix(prepared[["omega_i"]])),
                   rep(as.integer((k - 1) * tt), 2))
})

test_that("time varying covariance coefficients widen the regressor matrix", {
  skip_without_matrix()

  k <- 3
  tt <- 4
  u <- matrix(1:(k * tt))
  omega_i <- Matrix::Matrix(diag(1:k, k))

  constant <- covar_prepare_data(u, omega_i, k, tt, FALSE)
  varying <- covar_prepare_data(u, omega_i, k, tt, TRUE)

  expect_equal(ncol(as.matrix(varying[["z"]])),
               ncol(as.matrix(constant[["z"]])) * tt)
  expect_equal(as.matrix(varying[["y"]]), as.matrix(constant[["y"]]))
})

test_that("coint_prepare_sur_data builds the cointegration regressors", {
  skip_without_matrix()

  k <- 2
  r <- 1
  # w holds one row per observation and one column per cointegration regressor.
  w <- matrix(1:6, 3, 2)
  alpha <- matrix(c(0.5, -0.3), k, r)

  prepared <- coint_prepare_sur_data(w, alpha, k, r, FALSE, FALSE)

  expect_named(prepared, c("z", "alpha"))
  expect_equal(prepared[["alpha"]], alpha)

  # Each period contributes a k x (r * k_beta) block alpha %x% w_t'.
  expected <- do.call(rbind, lapply(seq_len(nrow(w)), function(i) {
    kronecker(alpha, matrix(w[i, ], 1))
  }))
  expect_equal(as.matrix(prepared[["z"]]), expected)
})
