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
