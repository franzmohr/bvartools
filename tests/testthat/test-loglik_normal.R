test_that("loglik_normal matches the multivariate normal density", {
  u <- matrix(c(1, 2, 3, 4, -1, 0.5), 2, 3)
  sigma <- matrix(c(2, 0.5, 0.5, 1), 2)

  expected <- vapply(seq_len(ncol(u)), function(i) {
    -0.5 * (2 * log(2 * pi) + log(det(sigma)) +
              t(u[, i]) %*% solve(sigma) %*% u[, i])
  }, numeric(1))

  expect_equal(as.numeric(loglik_normal(u, sigma)), expected)
})

test_that("loglik_normal returns one value per observation", {
  u <- matrix(stats::rnorm(30), 3, 10)

  expect_identical(dim(loglik_normal(u, diag(3))), c(10L, 1L))
})

test_that("a standard normal density is evaluated correctly", {
  u <- matrix(c(0, 0), 2, 1)

  # The density of N(0, I_2) at the origin.
  expect_equal(as.numeric(loglik_normal(u, diag(2))), -log(2 * pi))
})

test_that("larger residuals are less likely", {
  small <- loglik_normal(matrix(c(0.1, 0.1)), diag(2))
  large <- loglik_normal(matrix(c(5, 5)), diag(2))

  expect_lt(as.numeric(large), as.numeric(small))
})

test_that("time varying error variances are accepted", {
  u <- matrix(c(1, 2, 3, 4), 2, 2)
  # A 2T x 2 stack of per-period variance matrices.
  sigma <- rbind(diag(2), diag(2, 2))
  result <- loglik_normal(u, sigma)

  expect_identical(dim(result), c(2L, 1L))
  expect_equal(as.numeric(result[1]),
               as.numeric(loglik_normal(matrix(u[, 1]), diag(2))))
  expect_equal(as.numeric(result[2]),
               as.numeric(loglik_normal(matrix(u[, 2]), diag(2, 2))))
})
