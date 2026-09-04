test_that("ssvs returns inclusion indicators and a precision matrix", {
  set.seed(1)
  a <- matrix(c(0.01, 2, -3, 0.001))
  m <- nrow(a)

  result <- ssvs(a, matrix(0.05, m), matrix(10, m), matrix(0.5, m))

  expect_named(result, c("v_i", "lambda"))
  expect_identical(dim(result[["v_i"]]), c(m, m))
  expect_identical(dim(result[["lambda"]]), c(m, 1L))
  expect_true(all(result[["lambda"]] %in% c(0, 1)))
})

test_that("ssvs assigns the spike or the slab precision", {
  set.seed(2)
  a <- matrix(c(0.01, 2))
  tau0 <- 0.05
  tau1 <- 10

  result <- ssvs(a, matrix(tau0, 2), matrix(tau1, 2), matrix(0.5, 2))
  expected <- ifelse(result[["lambda"]] == 1, 1 / tau1^2, 1 / tau0^2)

  expect_equal(diag(result[["v_i"]]), as.numeric(expected))
  # Off-diagonal elements stay at zero.
  expect_true(all(result[["v_i"]][upper.tri(result[["v_i"]])] == 0))
})

test_that("ssvs excludes coefficients near zero and keeps large ones", {
  set.seed(3)
  # A spike this tight makes the classification effectively deterministic.
  a <- matrix(c(0, 5))
  result <- ssvs(a, matrix(1e-8, 2), matrix(10, 2), matrix(0.5, 2))

  expect_equal(as.numeric(result[["lambda"]]), c(0, 1))
})

test_that("ssvs can restrict selection to a subset of coefficients", {
  set.seed(4)
  a <- matrix(c(0, 0, 5, 5))
  result <- ssvs(a, matrix(1e-8, 4), matrix(10, 4), matrix(0.5, 4),
                 include = matrix(1:2))

  # Coefficients outside `include` are always kept in the model.
  expect_equal(as.numeric(result[["lambda"]]), c(0, 0, 1, 1))
})

test_that("post_bvs returns the inclusion matrix", {
  set.seed(5)

  k <- 2
  tt <- 30
  m <- k * k
  x <- matrix(stats::rnorm(k * tt), k, tt)
  # z is taken dense while lambda and sigma_i have to be sparse.
  z <- kronecker(t(x), diag(k))
  a <- matrix(c(0.5, 0, 0, 0.5))
  y <- matrix(z %*% a + stats::rnorm(k * tt))
  lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)
  sigma_i <- Matrix::Matrix(kronecker(diag(tt), diag(k)), sparse = TRUE)

  result <- as.matrix(post_bvs(y, z, a, k, m, lambda, sigma_i, matrix(0.5, m)))

  expect_equal(dim(result), c(m, m))
  expect_true(all(diag(result) %in% c(0, 1)))
  # Everything off the diagonal stays zero.
  expect_true(all(result[upper.tri(result)] == 0))
  expect_true(all(result[lower.tri(result)] == 0))
})

test_that("post_bvs keeps coefficients that clearly matter", {
  set.seed(6)

  k <- 2
  tt <- 400
  m <- k * k
  x <- matrix(stats::rnorm(k * tt), k, tt)
  z <- kronecker(t(x), diag(k))
  a <- matrix(c(3, 0, 0, 3))
  y <- matrix(z %*% a + stats::rnorm(k * tt) * 0.1)
  lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)
  sigma_i <- Matrix::Matrix(kronecker(diag(tt), diag(1 / 0.01, k)),
                            sparse = TRUE)

  result <- as.matrix(post_bvs(y, z, a, k, m, lambda, sigma_i, matrix(0.5, m)))

  # With a strong signal the two non-zero coefficients stay in the model.
  expect_equal(diag(result)[c(1, 4)], c(1, 1))
})

test_that("post_bvs can restrict selection to a subset of coefficients", {
  set.seed(7)

  k <- 2
  tt <- 50
  m <- k * k
  x <- matrix(stats::rnorm(k * tt), k, tt)
  z <- kronecker(t(x), diag(k))
  a <- matrix(c(0, 0, 0, 0))
  y <- matrix(stats::rnorm(k * tt))
  lambda <- Matrix::Matrix(diag(1, m), sparse = TRUE)
  sigma_i <- Matrix::Matrix(kronecker(diag(tt), diag(k)), sparse = TRUE)

  result <- as.matrix(post_bvs(y, z, a, k, m, lambda, sigma_i, matrix(0.5, m),
                               include = matrix(1:2)))

  # Positions outside `include` keep the value they came in with.
  expect_equal(diag(result)[3:4], c(1, 1))
})
