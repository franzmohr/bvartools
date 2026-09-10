test_that("fevd returns shares for every variable and horizon", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)
  spec <- fx_var_fitted()[["model"]]

  expect_s3_class(decomposition, "bvarfevd")
  expect_s3_class(decomposition, "ts")
  expect_identical(nrow(decomposition), 6L)
  expect_identical(colnames(decomposition), spec[["endogen"]])
  expect_equal(stats::start(decomposition)[1], 0)
})

test_that("the decomposition is a set of shares that add up", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)

  expect_true(all(decomposition >= 0))
  expect_true(all(decomposition <= 1))
  expect_equal(as.numeric(rowSums(decomposition)), rep(1, nrow(decomposition)))
})

test_that("an orthogonalised decomposition starts from the Cholesky order", {
  decomposition <- fevd(fx_var_fitted(), response = "invest", n_ahead = 3)

  # invest is ordered first, so on impact only its own shock explains it.
  expect_equal(as.numeric(decomposition[1, ]), c(1, 0, 0))
})

test_that("the generalised decomposition is also a set of shares", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 4,
                        type = "gir", normalise_gir = TRUE)

  expect_true(all(decomposition >= 0))
  expect_equal(as.numeric(rowSums(decomposition)), rep(1, nrow(decomposition)))
})

test_that("an unnormalised generalised decomposition need not add up to one", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 4,
                        type = "gir", normalise_gir = FALSE)

  expect_true(all(decomposition >= 0))
  expect_false(isTRUE(all.equal(as.numeric(rowSums(decomposition)),
                                rep(1, nrow(decomposition)))))
})

test_that("a structural model rejects the reduced-form types", {
  for (type in c("oir", "gir")) {
    expect_error(fevd(fx_svar_fitted(), response = "cons", n_ahead = 3,
                      type = type),
                 "not defined for a structural model")
  }
})

test_that("the structural decomposition accounts for the structural variances", {
  # A recursive structural model and the reduced form it implies. The two
  # describe the same process, so the structural decomposition of the one has to
  # equal the orthogonalised decomposition of the other -- which only holds if
  # the structural variances enter the decomposition.
  k <- 3
  A0 <- diag(k)
  A0[lower.tri(A0)] <- c(0.4, -0.7, 1.3)
  Sigma <- diag(c(0.09, 0.64, 2.25)) # deliberately far from the identity
  A <- matrix(c(0.5, 0.1, -0.2,
                0.0, 0.3, 0.4,
                0.1, -0.1, 0.6), k, k)

  A0_inv <- solve(A0)
  reduced <- list(A = A0_inv %*% A, Sigma = A0_inv %*% Sigma %*% t(A0_inv))

  for (response in 1:k) {
    structural <- bvartools:::.vardecomp(
      list(A = A0_inv %*% A, Sigma = Sigma, A0 = A0),
      h = 5, type = "sir", response = response)

    expect_equal(structural, bvartools:::.vardecomp(reduced, h = 5, type = "oir",
                                                    response = response))
    expect_equal(rowSums(structural), rep(1, nrow(structural)))
  }
})

test_that("a structural decomposition of a fitted model adds up", {
  decomposition <- fevd(fx_svar_fitted(), response = "cons", n_ahead = 4,
                        type = "sir")

  expect_s3_class(decomposition, "bvarfevd")
  expect_true(all(decomposition >= 0))
  expect_equal(as.numeric(rowSums(decomposition)), rep(1, nrow(decomposition)))
})

test_that("unknown response variables are rejected", {
  expect_error(fevd(fx_var_fitted(), response = "nonexistent"))
})

test_that("a converted VEC model can be decomposed", {
  decomposition <- fevd(vec_to_var(fx_vec_fitted()), response = "Dp",
                        n_ahead = 4)

  expect_s3_class(decomposition, "bvarfevd")
  expect_equal(as.numeric(rowSums(decomposition)), rep(1, nrow(decomposition)))
})

test_that("all variables are shown when no maximum is given", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 3)
  spec <- fx_var_fitted()[["model"]]

  expect_identical(colnames(decomposition), spec[["endogen"]])
})

test_that("max_groups keeps the largest contributions and pools the rest", {
  full <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)
  limited <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5,
                  max_groups = 2)

  expect_identical(ncol(limited), 2L)

  # The variable with the largest total contribution is kept by name, the
  # others end up in "Other".
  largest <- colnames(full)[which.max(colSums(full))]
  expect_identical(colnames(limited), c(largest, "Other"))
  expect_equal(as.numeric(limited[, largest]), as.numeric(full[, largest]))
  expect_equal(as.numeric(limited[, "Other"]),
               as.numeric(rowSums(full[, colnames(full) != largest])))
})

test_that("pooling leaves the row sums of the decomposition untouched", {
  limited <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5,
                  max_groups = 2)

  expect_equal(as.numeric(rowSums(limited)), rep(1, nrow(limited)))
  expect_s3_class(limited, "bvarfevd")
})

test_that("a maximum of at least the number of variables changes nothing", {
  full <- fevd(fx_var_fitted(), response = "cons", n_ahead = 3)
  limited <- fevd(fx_var_fitted(), response = "cons", n_ahead = 3,
                  max_groups = ncol(full) + 1)

  expect_equal(limited, full)
})

test_that("an implausible maximum is rejected", {
  expect_error(fevd(fx_var_fitted(), response = "cons", max_groups = 0),
               "positive integer")
  expect_error(fevd(fx_var_fitted(), response = "cons", max_groups = c(1, 2)),
               "positive integer")
  expect_error(fevd(fx_var_fitted(), response = "cons", max_groups = "two"),
               "positive integer")
})
