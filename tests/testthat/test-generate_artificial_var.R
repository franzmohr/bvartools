test_that("the generated series has the requested shape", {
  set.seed(1)
  artificial <- generate_artificial_var(nobs = 60, k = 3, p = 2)

  expect_named(artificial, c("data", "params"))
  expect_s3_class(artificial[["data"]], "ts")
  expect_identical(dim(artificial[["data"]]), c(60L, 3L))
  expect_true(all(is.finite(artificial[["data"]])))
})

test_that("the coefficient matrix has one block per lag", {
  set.seed(2)
  artificial <- generate_artificial_var(nobs = 50, k = 2, p = 3)

  expect_identical(dim(artificial[["params"]][["a_coef"]]), c(2L, 6L))
  expect_identical(dim(artificial[["params"]][["u_sigma"]]), c(2L, 2L))
})

test_that("a_zeros controls how many coefficients are zero", {
  set.seed(3)
  none <- generate_artificial_var(nobs = 40, k = 3, p = 2, a_zeros = 0)
  set.seed(3)
  all_zero <- generate_artificial_var(nobs = 40, k = 3, p = 2, a_zeros = 1)

  expect_true(all(none[["params"]][["a_coef"]] != 0))
  expect_true(all(all_zero[["params"]][["a_coef"]] == 0))
})

test_that("coefficients stay inside the requested range", {
  set.seed(4)
  artificial <- generate_artificial_var(nobs = 40, k = 3, p = 2, a_zeros = 0,
                                        a_range = c(-0.2, 0.2))
  coefficients <- artificial[["params"]][["a_coef"]]

  expect_true(all(abs(coefficients) <= 0.2))
})

test_that("a constant is appended to the coefficient matrix", {
  set.seed(5)
  artificial <- generate_artificial_var(nobs = 40, k = 2, p = 1, const = TRUE,
                                        range_const = c(1, 2))
  coefficients <- artificial[["params"]][["a_coef"]]

  # k * p lag coefficients plus one intercept column.
  expect_identical(dim(coefficients), c(2L, 3L))
  intercept <- coefficients[, 3]
  expect_true(all(intercept >= 1 & intercept <= 2))
})

test_that("a presample is discarded from the output", {
  set.seed(6)
  artificial <- generate_artificial_var(nobs = 40, k = 2, p = 1,
                                        presample = 20)

  expect_identical(nrow(artificial[["data"]]), 40L)
})

test_that("the generated data can be estimated", {
  set.seed(7)
  artificial <- generate_artificial_var(nobs = 80, k = 2, p = 1, a_zeros = 0)

  model <- create_bvarmodel(artificial[["data"]], p = 1,
                            deterministic = "none",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  model <- add_initial_values(model)
  set.seed(8)
  model <- add_posterior_coefficients(model)

  expect_identical(dim(model[["posterior"]][["a"]][["coeffs"]]), c(10L, 4L))
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
})

test_that("invalid arguments are rejected", {
  expect_error(generate_artificial_var(a_zeros = 1.5), "between 0 and 1")
  expect_error(generate_artificial_var(a_zeros = -1), "between 0 and 1")
  expect_error(generate_artificial_var(const = TRUE),
               "'range_const' must be specified")
})
