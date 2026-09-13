test_that("the generated series has the requested shape", {
  set.seed(1)
  artificial <- generate_artificial_var(nobs = 60, k = 3, p = 2)

  expect_named(artificial, c("data", "params"))
  expect_s3_class(artificial[["data"]], "ts")
  expect_identical(dim(artificial[["data"]]), c(60L, 3L))
  expect_identical(colnames(artificial[["data"]]), c("var1", "var2", "var3"))
  expect_true(all(is.finite(artificial[["data"]])))
  expect_named(artificial[["params"]], c("a_coef", "psi_coef", "u_omega", "u_sigma"))
})

test_that("the coefficient matrix has one block per lag", {
  set.seed(2)
  artificial <- generate_artificial_var(nobs = 50, k = 2, p = 3)

  expect_identical(dim(artificial[["params"]][["a_coef"]]), c(2L, 6L))
  expect_identical(colnames(artificial[["params"]][["a_coef"]]),
                   c("var1.1", "var2.1", "var1.2", "var2.2", "var1.3", "var2.3"))
  expect_identical(dim(artificial[["params"]][["u_sigma"]]), c(2L, 2L))

  # Every coefficient is drawn on its own rather than recycled over the lags.
  set.seed(2)
  dense <- generate_artificial_var(nobs = 50, k = 3, p = 2, a_zeros = 0)
  expect_length(unique(as.vector(dense[["params"]][["a_coef"]])), 18L)
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
                                        range_a = c(0.2, -0.2))
  coefficients <- artificial[["params"]][["a_coef"]]

  expect_true(all(abs(coefficients) <= 0.2))
})

test_that("drawn coefficients describe a stable process", {
  set.seed(5)
  for (p in 1:4) {
    artificial <- generate_artificial_var(nobs = 20, k = 3, p = p, a_zeros = 0,
                                          range_a = c(-0.6, 0.6))
    expect_true(.artificial_var_stable(artificial[["params"]][["a_coef"]], 3, p))
  }
})

test_that("deterministic terms are appended to the coefficient matrix", {
  set.seed(6)
  artificial <- generate_artificial_var(nobs = 40, k = 2, p = 1, deterministic = "both",
                                        range_const = c(1, 2), range_trend = c(3, 4))
  coefficients <- artificial[["params"]][["a_coef"]]

  # k * p lag coefficients plus an intercept and a trend column.
  expect_identical(dim(coefficients), c(2L, 4L))
  expect_identical(colnames(coefficients), c("var1.1", "var2.1", "const", "trend"))
  expect_true(all(coefficients[, "const"] >= 1 & coefficients[, "const"] <= 2))
  expect_true(all(coefficients[, "trend"] >= 3 & coefficients[, "trend"] <= 4))

  set.seed(6)
  none <- generate_artificial_var(nobs = 40, k = 2, p = 0)
  expect_true("a_coef" %in% names(none[["params"]]))
  expect_null(none[["params"]][["a_coef"]])
})

test_that("the trend is aligned with the trend of create_bvarmodel", {
  set.seed(7)
  artificial <- generate_artificial_var(nobs = 30, k = 2, p = 2, deterministic = "trend",
                                        a_zeros = 1, range_trend = c(1, 1),
                                        range_variance = c(1e-12, 1e-12))
  model <- create_bvarmodel(artificial[["data"]], p = 2, deterministic = "trend")
  x <- model[["data"]][["train"]][["x"]]

  # Without lags or noise the series is its trend.
  expect_equal(as.numeric(artificial[["data"]][3:30, 1]), as.numeric(x[, "trend"]),
               tolerance = 1e-5)
})

test_that("a presample is discarded from the output", {
  set.seed(8)
  artificial <- generate_artificial_var(nobs = 40, k = 2, p = 1,
                                        presample = 20)

  expect_identical(nrow(artificial[["data"]]), 40L)
})

test_that("the error covariance follows from psi and omega", {
  set.seed(9)
  artificial <- generate_artificial_var(nobs = 40, k = 3, p = 1,
                                        range_variance = c(0.5, 2), range_psi = c(-1, 1))
  params <- artificial[["params"]]
  psi_inv <- solve(params[["psi_coef"]])

  expect_true(all(params[["psi_coef"]][upper.tri(params[["psi_coef"]])] == 0))
  expect_equal(params[["u_sigma"]], psi_inv %*% params[["u_omega"]] %*% t(psi_inv))
})

test_that("time varying parameters are returned as paths", {
  set.seed(10)
  artificial <- generate_artificial_var(nobs = 50, k = 2, p = 2, deterministic = "const",
                                        tvp = TRUE, a_zeros = 0.3, range_psi = c(0.5, 0.5),
                                        range_variance_state = c(0.001, 0.01))
  params <- artificial[["params"]]
  a_path <- params[["a_coef"]]

  expect_identical(dim(a_path), c(2L, 5L, 50L))
  expect_identical(dim(params[["psi_coef"]]), c(2L, 2L, 50L))
  expect_identical(dim(params[["u_sigma"]]), c(2L, 2L, 50L))
  expect_identical(dim(params[["u_omega"]]), c(2L, 2L))
  expect_identical(dim(params[["a_state_variance"]]), c(2L, 5L))

  # Non-zero coefficients move, zero coefficients and their state variances stay zero.
  zeros <- a_path[, , 1] == 0
  expect_true(any(zeros))
  expect_true(all(a_path[, , 50][zeros] == 0))
  expect_true(all(params[["a_state_variance"]][zeros] == 0))
  expect_true(all(apply(a_path, c(1, 2), stats::sd)[!zeros] > 0))
  expect_gt(stats::sd(params[["psi_coef"]][2, 1, ]), 0)
  expect_true(all(params[["psi_coef"]][1, 2, ] == 0))

  # The process is stable in every period.
  for (t in 1:50) {
    expect_true(.artificial_var_stable(a_path[, 1:4, t], 2, 2))
  }
})

test_that("stochastic volatility is returned as paths", {
  set.seed(11)
  artificial <- generate_artificial_var(nobs = 50, k = 3, p = 1, sv = TRUE)
  params <- artificial[["params"]]
  omega <- params[["u_omega"]]

  expect_identical(dim(params[["a_coef"]]), c(3L, 3L))
  expect_identical(dim(omega), c(3L, 3L, 50L))
  expect_identical(dim(params[["u_sigma"]]), c(3L, 3L, 50L))
  expect_length(params[["u_state_variance"]], 3)
  expect_true(all(omega[1, 2, ] == 0))
  expect_true(all(apply(omega, 1:2, stats::sd)[diag(3) == 1] > 0))
  expect_null(params[["a_state_variance"]])
})

test_that("the data are generated from the returned parameters", {
  set.seed(12)
  nobs <- 3000
  p <- 2
  artificial <- generate_artificial_var(nobs = nobs, k = 2, p = p, deterministic = "both",
                                        tvp = TRUE, sv = TRUE, range_psi = c(-0.5, 0.5),
                                        range_variance = c(0.5, 2))
  y <- artificial[["data"]]
  params <- artificial[["params"]]

  # Standardised structural errors of the returned periods after the first p.
  errors <- sapply((p + 1):nobs, function(t) {
    x <- c(y[t - 1, ], y[t - 2, ], 1, t - p)
    u <- y[t, ] - params[["a_coef"]][, , t] %*% x
    (params[["psi_coef"]][, , t] %*% u) / sqrt(diag(params[["u_omega"]][, , t]))
  })

  expect_equal(rowMeans(errors), c(0, 0), tolerance = 0.1)
  expect_equal(apply(errors, 1, stats::var), c(1, 1), tolerance = 0.1)
  expect_lt(abs(stats::cor(errors[1, ], errors[2, ])), 0.1)
})

test_that("structural models are generated from A0", {
  set.seed(13)
  nobs <- 3000
  artificial <- generate_artificial_var(nobs = nobs, k = 3, p = 1, deterministic = "const",
                                        structural = TRUE, tvp = TRUE, sv = TRUE,
                                        range_a0 = c(-1, 1), range_variance = c(0.5, 2))
  y <- artificial[["data"]]
  params <- artificial[["params"]]
  a0 <- params[["a0_coef"]]

  expect_named(params, c("a_coef", "a0_coef", "psi_coef", "u_omega", "u_sigma",
                         "a_state_variance", "a0_state_variance", "u_state_variance"))
  expect_null(params[["psi_coef"]])
  expect_identical(dim(a0), c(3L, 3L, as.integer(nobs)))
  expect_true(all(a0[2, 2, ] == 1 & a0[1, 3, ] == 0))
  expect_true(all(params[["a0_state_variance"]][upper.tri(diag(3), diag = TRUE)] == 0))
  expect_gt(stats::sd(a0[3, 2, ]), 0)
  expect_equal(params[["u_sigma"]], params[["u_omega"]])

  errors <- sapply(2:nobs, function(t) {
    u <- a0[, , t] %*% y[t, ] - params[["a_coef"]][, , t] %*% c(y[t - 1, ], 1)
    u / sqrt(diag(params[["u_omega"]][, , t]))
  })

  expect_equal(apply(errors, 1, stats::var), rep(1, 3), tolerance = 0.1)
  expect_lt(max(abs(stats::cor(t(errors))[lower.tri(diag(3))])), 0.1)

  # The reduced form is stable in every period.
  for (t in c(1, nobs)) {
    reduced <- solve(a0[, , t], params[["a_coef"]][, 1:3, t])
    expect_true(.artificial_var_stable(reduced, 3, 1))
  }
})

test_that("the generated data can be estimated", {
  set.seed(14)
  artificial <- generate_artificial_var(nobs = 80, k = 2, p = 1, a_zeros = 0)

  model <- create_bvarmodel(artificial[["data"]], p = 1,
                            deterministic = "none",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  model <- add_initial_values(model)
  set.seed(15)
  model <- add_posterior_coefficients(model)

  expect_identical(dim(model[["posterior"]][["a"]][["coeffs"]]), c(10L, 4L))
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
})

test_that("a univariate model has no psi", {
  set.seed(16)
  artificial <- generate_artificial_var(nobs = 30, k = 1, p = 1, tvp = TRUE, sv = TRUE)

  expect_true("psi_coef" %in% names(artificial[["params"]]))
  expect_null(artificial[["params"]][["psi_coef"]])
  expect_null(artificial[["params"]][["psi_state_variance"]])
  expect_identical(dim(artificial[["params"]][["u_sigma"]]), c(1L, 1L, 30L))
})

test_that("stable = FALSE allows integrated series", {
  set.seed(17)
  artificial <- generate_artificial_var(nobs = 30, k = 1, p = 1, range_a = c(1, 1),
                                        a_zeros = 0, stable = FALSE, presample = 0)

  expect_equal(unname(artificial[["params"]][["a_coef"]]), matrix(1))
  expect_error(generate_artificial_var(k = 1, p = 1, range_a = c(1, 1), a_zeros = 0),
               "No stable coefficient matrix")
})

test_that("invalid arguments are rejected", {
  expect_error(generate_artificial_var(a_zeros = 1.5), "between 0 and 1")
  expect_error(generate_artificial_var(a_zeros = -1), "between 0 and 1")
  expect_error(generate_artificial_var(deterministic = "season"), "'deterministic' must be one of")
  expect_error(generate_artificial_var(p = -1), "'p' must be a whole number")
  expect_error(generate_artificial_var(nobs = 2.5), "'nobs' must be a whole number")
  expect_error(generate_artificial_var(tvp = NA), "'tvp' must be TRUE or FALSE")
  expect_error(generate_artificial_var(range_variance = c(0, 1)), "must be positive")
  expect_error(generate_artificial_var(range_variance_sv = c(-1, 1)), "must not be negative")
  expect_error(generate_artificial_var(range_psi = 1), "two finite elements")
  expect_error(generate_artificial_var(structural = TRUE, range_psi = c(-1, 1)), "structural models")
})
