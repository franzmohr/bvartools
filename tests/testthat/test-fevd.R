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
