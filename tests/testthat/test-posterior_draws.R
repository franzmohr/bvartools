# The samplers are random, so they are checked against the analytic posterior
# they are supposed to draw from: either by making the posterior degenerate, or
# by averaging enough draws.

test_that("post_normal returns one value per coefficient", {
  set.seed(1)
  y <- matrix(stats::rnorm(20), 2, 10)
  x <- matrix(stats::rnorm(30), 3, 10)

  draw <- post_normal(y, x, diag(2), matrix(0, 6), diag(6))
  expect_identical(dim(draw), c(6L, 1L))
  expect_true(all(is.finite(draw)))
})

test_that("post_normal collapses to the prior mean under a tight prior", {
  set.seed(2)
  y <- matrix(stats::rnorm(20), 2, 10)
  x <- matrix(stats::rnorm(20), 2, 10)
  prior_mean <- matrix(c(0.5, -1, 2, 0.25))

  draw <- post_normal(y, x, diag(2), prior_mean, diag(1e12, 4))
  expect_equal(as.numeric(draw), as.numeric(prior_mean), tolerance = 1e-4)
})

test_that("post_normal draws around the analytic posterior mean", {
  set.seed(3)
  k <- 2
  tt <- 60
  x <- matrix(stats::rnorm(k * tt), k, tt)
  y <- matrix(c(0.4, -0.2, 0.1, 0.6), k) %*% x +
    matrix(stats::rnorm(k * tt), k, tt)

  sigma_i <- diag(k)
  prior_mean <- matrix(0, k * k)
  prior_precision <- diag(0.5, k * k)

  posterior_precision <- prior_precision + kronecker(tcrossprod(x), sigma_i)
  posterior_mean <- solve(posterior_precision) %*%
    (prior_precision %*% prior_mean + matrix(sigma_i %*% tcrossprod(y, x)))

  draws <- replicate(4000,
                     post_normal(y, x, sigma_i, prior_mean, prior_precision))
  expect_equal(rowMeans(draws[, 1, ]), as.numeric(posterior_mean),
               tolerance = 0.05)
  expect_equal(stats::var(t(draws[, 1, ])), solve(posterior_precision),
               tolerance = 0.1)
})

test_that("post_normal is reproducible from the seed", {
  y <- matrix(stats::rnorm(20), 2, 10)
  x <- matrix(stats::rnorm(20), 2, 10)

  set.seed(4)
  first <- post_normal(y, x, diag(2), matrix(0, 4), diag(4))
  set.seed(4)
  second <- post_normal(y, x, diag(2), matrix(0, 4), diag(4))
  set.seed(5)
  third <- post_normal(y, x, diag(2), matrix(0, 4), diag(4))

  expect_identical(first, second)
  expect_false(isTRUE(all.equal(first, third)))
})

test_that("post_normal_sur agrees with post_normal on the same model", {
  set.seed(6)
  k <- 2
  tt <- 40
  x <- matrix(stats::rnorm(k * tt), k, tt)
  y <- matrix(stats::rnorm(k * tt), k, tt)
  z <- kronecker(t(x), diag(k))

  prior_mean <- matrix(c(1, -1, 0.5, 0))
  # A degenerate prior makes both samplers return the same point.
  expect_equal(
    as.numeric(post_normal(y, x, diag(k), prior_mean, diag(1e12, k * k))),
    as.numeric(post_normal_sur(y, z, diag(k), prior_mean, diag(1e12, k * k))),
    tolerance = 1e-4
  )
})

test_that("post_normal_sur can use the singular value decomposition", {
  set.seed(7)
  k <- 2
  tt <- 30
  z <- kronecker(matrix(stats::rnorm(k * tt), tt, k), diag(k))
  y <- matrix(stats::rnorm(k * tt), k, tt)
  prior_mean <- matrix(c(0.3, -0.4, 0.1, 0.2))

  expect_equal(
    as.numeric(post_normal_sur(y, z, diag(k), prior_mean, diag(1e12, k * k),
                               svd = TRUE)),
    as.numeric(prior_mean),
    tolerance = 1e-4
  )
})

test_that("post_gamma_measurement_variance draws positive variances", {
  set.seed(8)
  k <- 4
  u <- matrix(stats::rnorm(k * 200))

  # The draw comes back as a sparse matrix.
  draw <- as.matrix(post_gamma_measurement_variance(
    u, matrix(1, k), matrix(0.0001, k), inverse = FALSE))

  expect_equal(dim(draw), c(k, k))
  expect_true(all(diag(draw) > 0))
  # It is a diagonal matrix.
  expect_true(all(draw[upper.tri(draw)] == 0))
})

test_that("the inverse option returns the precision", {
  set.seed(9)
  k <- 3
  u <- matrix(stats::rnorm(k * 500))
  arguments <- list(u = u, shape_prior = matrix(1, k),
                    rate_prior = matrix(0.0001, k))

  set.seed(10)
  variance <- as.matrix(do.call(post_gamma_measurement_variance,
                                c(arguments, list(inverse = FALSE))))
  set.seed(10)
  precision <- as.matrix(do.call(post_gamma_measurement_variance,
                                 c(arguments, list(inverse = TRUE))))

  expect_equal(diag(precision), 1 / diag(variance))
})

test_that("post_gamma_measurement_variance recovers a known variance", {
  set.seed(11)
  k <- 2
  tt <- 5000
  u <- matrix(stats::rnorm(k * tt, sd = c(2, 0.5)))

  draw <- as.matrix(post_gamma_measurement_variance(
    u, matrix(0.01, k), matrix(0.01, k), inverse = FALSE))
  expect_equal(diag(draw), c(4, 0.25), tolerance = 0.15)
})

test_that("post_gamma_state_variance recovers a known state variance", {
  set.seed(12)
  k <- 2
  tt <- 5000
  state_sd <- c(0.2, 0.05)

  states <- matrix(0, k, tt + 1)
  for (i in 2:(tt + 1)) {
    states[, i] <- states[, i - 1] + stats::rnorm(k, 0, state_sd)
  }
  initial <- matrix(states[, 1])
  path <- matrix(states[, -1])

  draw <- as.matrix(post_gamma_state_variance(
    path, initial, matrix(0.01, k), matrix(0.01, k), inverse = FALSE))

  expect_equal(dim(draw), c(k, k))
  expect_equal(diag(draw), state_sd^2, tolerance = 0.15)
})
