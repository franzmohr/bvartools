# Number of roots of the VAR representation of a VEC model on the unit circle
unit_roots <- function(pi, gamma, k, p) {
  a <- matrix(0, k, k * p)
  a[, 1:k] <- diag(1, k) + pi
  if (p > 1) {
    for (i in 1:(p - 1)) {
      a[, (i - 1) * k + 1:k] <- a[, (i - 1) * k + 1:k] + gamma[, (i - 1) * k + 1:k]
      a[, i * k + 1:k] <- a[, i * k + 1:k] - gamma[, (i - 1) * k + 1:k]
    }
  }
  roots <- Mod(eigen(vars_companion(a, k, p), only.values = TRUE)$values)
  list(unit = sum(abs(roots - 1) < 1e-6), max_other = max(c(0, roots[abs(roots - 1) >= 1e-6])))
}

vars_companion <- function(a, k, p) {
  companion <- matrix(0, k * p, k * p)
  companion[1:k, ] <- a
  if (p > 1) {
    companion[k + 1:(k * (p - 1)), 1:(k * (p - 1))] <- diag(1, k * (p - 1))
  }
  companion
}

test_that("the generated series has the requested shape", {
  set.seed(1)
  artificial <- generate_artificial_vec(nobs = 60, k = 3, p = 2, r = 1)
  params <- artificial[["params"]]

  expect_named(artificial, c("data", "params"))
  expect_s3_class(artificial[["data"]], "ts")
  expect_identical(dim(artificial[["data"]]), c(60L, 3L))
  expect_true(all(is.finite(artificial[["data"]])))
  expect_named(params, c("alpha", "beta", "pi", "gamma", "c", "psi_coef", "u_omega", "u_sigma"))
  expect_identical(dim(params[["alpha"]]), c(3L, 1L))
  expect_identical(dim(params[["beta"]]), c(3L, 1L))
  expect_identical(unname(params[["beta"]][1, 1]), 1)
  expect_equal(params[["pi"]], params[["alpha"]] %*% t(params[["beta"]]), ignore_attr = TRUE)
  expect_identical(colnames(params[["pi"]]), c("l.var1", "l.var2", "l.var3"))
  expect_identical(colnames(params[["gamma"]]), c("d.var1.l01", "d.var2.l01", "d.var3.l01"))
  expect_null(params[["c"]])
})

test_that("the process has k - r unit roots", {
  set.seed(2)
  specs <- list(c(k = 3, p = 1, r = 1), c(k = 3, p = 2, r = 2), c(k = 2, p = 3, r = 1),
                c(k = 3, p = 2, r = 0), c(k = 2, p = 2, r = 2))
  for (spec in specs) {
    k <- spec[["k"]]
    p <- spec[["p"]]
    r <- spec[["r"]]
    params <- generate_artificial_vec(nobs = 20, k = k, p = p, r = r, gamma_zeros = 0)[["params"]]
    pi <- if (r > 0) params[["pi"]] else matrix(0, k, k)
    roots <- unit_roots(pi, params[["gamma"]], k, p)

    expect_identical(roots[["unit"]], as.integer(k - r))
    expect_lt(roots[["max_other"]], 1)
    if (r == 0) {
      expect_null(params[["alpha"]])
      expect_null(params[["pi"]])
    }
  }
})

test_that("deterministic terms enter the cointegration term or the short-run part", {
  set.seed(3)
  artificial <- generate_artificial_vec(nobs = 40, k = 2, p = 1, r = 1,
                                        const = "restricted", trend = "unrestricted",
                                        range_const = c(1, 2), range_trend = c(3, 4))
  params <- artificial[["params"]]

  expect_identical(rownames(params[["beta"]]), c("l.var1", "l.var2", "const"))
  expect_true(params[["beta"]]["const", 1] >= 1 & params[["beta"]]["const", 1] <= 2)
  expect_identical(colnames(params[["c"]]), "trend")
  expect_true(all(params[["c"]] >= 3 & params[["c"]] <= 4))
  expect_null(params[["gamma"]])

  expect_error(generate_artificial_vec(r = 0, const = "restricted"), "at least 1")
})

test_that("the trend is aligned with the trend of create_bvecmodel", {
  set.seed(4)
  artificial <- generate_artificial_vec(nobs = 30, k = 2, p = 2, r = 0, trend = "unrestricted",
                                        gamma_zeros = 1, range_trend = c(1, 1),
                                        range_variance = c(1e-12, 1e-12))
  model <- create_bvecmodel(artificial[["data"]], p = 2, r = 0, trend = "unrestricted")
  x <- model[["data"]][["train"]][["x"]]

  # Without short-run dynamics or noise the differences are the trend.
  expect_equal(as.numeric(model[["data"]][["train"]][["y"]][, 1]), as.numeric(x[, ncol(x)]),
               tolerance = 1e-5)
})

test_that("the data are generated from the returned parameters", {
  set.seed(5)
  nobs <- 3000
  artificial <- generate_artificial_vec(nobs = nobs, k = 3, p = 2, r = 1,
                                        const = "restricted", trend = "unrestricted",
                                        tvp = TRUE, sv = TRUE, range_psi = c(-0.5, 0.5),
                                        range_variance = c(0.5, 2))
  y <- artificial[["data"]]
  params <- artificial[["params"]]

  errors <- sapply(3:nobs, function(t) {
    w <- c(y[t - 1, ], 1)
    u <- y[t, ] - y[t - 1, ] - params[["pi"]][, , t] %*% w -
      params[["gamma"]][, , t] %*% (y[t - 1, ] - y[t - 2, ]) -
      params[["c"]][, "trend", t] * (t - 2)
    (params[["psi_coef"]][, , t] %*% u) / sqrt(diag(params[["u_omega"]][, , t]))
  })

  expect_equal(rowMeans(errors), rep(0, 3), tolerance = 0.1)
  expect_equal(apply(errors, 1, stats::var), rep(1, 3), tolerance = 0.1)
  expect_lt(max(abs(stats::cor(t(errors))[lower.tri(diag(3))])), 0.1)
})

test_that("time varying parameters are returned as paths", {
  set.seed(6)
  nobs <- 50
  artificial <- generate_artificial_vec(nobs = nobs, k = 3, p = 2, r = 2, const = "unrestricted",
                                        tvp = TRUE, gamma_zeros = 0.5,
                                        range_variance_state = c(0.001, 0.01))
  params <- artificial[["params"]]

  expect_identical(dim(params[["alpha"]]), c(3L, 2L, 50L))
  expect_identical(dim(params[["beta"]]), c(3L, 2L, 50L))
  expect_identical(dim(params[["pi"]]), c(3L, 3L, 50L))
  expect_identical(dim(params[["c"]]), c(3L, 1L, 50L))
  expect_identical(dim(params[["beta_state_variance"]]), c(3L, 2L))

  # The normalisation of beta is fixed, the free elements move.
  expect_true(all(params[["beta"]][1:2, , ] == array(diag(1, 2), c(2, 2, nobs))))
  expect_true(all(params[["beta_state_variance"]][1:2, ] == 0))
  expect_gt(stats::sd(params[["beta"]][3, 1, ]), 0)

  zeros <- params[["gamma"]][, , 1] == 0
  expect_true(any(zeros))
  expect_true(all(params[["gamma"]][, , nobs][zeros] == 0))

  for (t in 1:nobs) {
    expect_equal(params[["pi"]][, , t], params[["alpha"]][, , t] %*% t(params[["beta"]][, , t]),
                 ignore_attr = TRUE)
    roots <- unit_roots(params[["pi"]][, , t], params[["gamma"]][, , t], 3, 2)
    expect_identical(roots[["unit"]], 1L)
    expect_lt(roots[["max_other"]], 1)
  }
})

test_that("structural models are generated from A0", {
  set.seed(7)
  nobs <- 3000
  artificial <- generate_artificial_vec(nobs = nobs, k = 3, p = 2, r = 1, const = "unrestricted",
                                        structural = TRUE, tvp = TRUE, sv = TRUE,
                                        range_a0 = c(-1, 1), range_variance = c(0.5, 2))
  y <- artificial[["data"]]
  params <- artificial[["params"]]
  a0 <- params[["a0_coef"]]

  expect_null(params[["psi_coef"]])
  expect_null(params[["psi_state_variance"]])
  expect_identical(dim(a0), c(3L, 3L, as.integer(nobs)))
  expect_true(all(a0[1, 1, ] == 1 & a0[1, 2, ] == 0 & a0[2, 3, ] == 0))
  expect_gt(stats::sd(a0[3, 1, ]), 0)
  expect_equal(params[["u_sigma"]], params[["u_omega"]])

  errors <- sapply(3:nobs, function(t) {
    u <- a0[, , t] %*% (y[t, ] - y[t - 1, ]) - params[["pi"]][, , t] %*% y[t - 1, ] -
      params[["gamma"]][, , t] %*% (y[t - 1, ] - y[t - 2, ]) - params[["c"]][, "const", t]
    u / sqrt(diag(params[["u_omega"]][, , t]))
  })

  expect_equal(apply(errors, 1, stats::var), rep(1, 3), tolerance = 0.1)
  expect_lt(max(abs(stats::cor(t(errors))[lower.tri(diag(3))])), 0.1)
})

test_that("the generated data can be estimated", {
  set.seed(8)
  artificial <- generate_artificial_vec(nobs = 100, k = 2, p = 2, r = 1, const = "unrestricted")

  model <- create_bvecmodel(artificial[["data"]], p = 2, r = 1, const = "unrestricted",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1))
  set.seed(9)
  model <- add_initial_values(model)
  model <- add_posterior_coefficients(model)

  expect_identical(nrow(model[["posterior"]][["beta"]][["coeffs"]]), 10L)
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
})

test_that("invalid arguments are rejected", {
  expect_error(generate_artificial_vec(r = 4), "must not be larger than 'k'")
  expect_error(generate_artificial_vec(p = 0), "'p' must be a whole number of at least 1")
  expect_error(generate_artificial_vec(const = "both"), "NULL, 'restricted' or 'unrestricted'")
  expect_error(generate_artificial_vec(gamma_zeros = 2), "between 0 and 1")
  expect_error(generate_artificial_vec(structural = TRUE, range_psi = c(-1, 1)), "structural models")
  expect_error(generate_artificial_vec(range_alpha = c(1, 1), range_beta = c(1, 1)),
               "No coefficients of a process")
})

test_that("the level shifts the series and the constant absorbs it", {
  level <- c(100, 50, 20)
  for (const in c("restricted", "unrestricted")) {
    set.seed(10)
    base <- generate_artificial_vec(nobs = 50, k = 3, p = 2, r = 1, const = const)
    set.seed(10)
    shifted <- generate_artificial_vec(nobs = 50, k = 3, p = 2, r = 1, const = const, level = level)
    y <- shifted[["data"]]
    params <- shifted[["params"]]

    expect_equal(unclass(y) - unclass(base[["data"]]), matrix(level, 50, 3, byrow = TRUE),
                 ignore_attr = TRUE)

    # The returned parameters generate the shifted series without errors beyond those of the base
    det <- if (const == "restricted") c(y[2, ], 1) else y[2, ]
    fit <- params[["pi"]] %*% det + params[["gamma"]] %*% (y[2, ] - y[1, ])
    if (const == "unrestricted") {
      fit <- fit + params[["c"]][, "const"]
    }
    det_base <- if (const == "restricted") c(base[["data"]][2, ], 1) else base[["data"]][2, ]
    fit_base <- base[["params"]][["pi"]] %*% det_base +
      base[["params"]][["gamma"]] %*% (base[["data"]][2, ] - base[["data"]][1, ])
    if (const == "unrestricted") {
      fit_base <- fit_base + base[["params"]][["c"]][, "const"]
    }
    expect_equal(fit, fit_base)
    expect_equal(params[["alpha"]], base[["params"]][["alpha"]])
  }

  # Time varying parameters are shifted in every period
  set.seed(11)
  base <- generate_artificial_vec(nobs = 20, k = 2, p = 1, r = 1, const = "restricted", tvp = TRUE,
                                  range_variance_state = c(0.001, 0.001))
  set.seed(11)
  shifted <- generate_artificial_vec(nobs = 20, k = 2, p = 1, r = 1, const = "restricted", tvp = TRUE,
                                     range_variance_state = c(0.001, 0.001), level = 10)
  b <- base[["params"]][["beta"]]
  expect_equal(shifted[["params"]][["beta"]]["const", 1, ], b["const", 1, ] - 10 * colSums(b[1:2, 1, ]))
  expect_equal(shifted[["params"]][["pi"]][, 1:2, ], base[["params"]][["pi"]][, 1:2, ])

  # Without cointegration no constant is needed
  set.seed(12)
  expect_equal(unname(generate_artificial_vec(nobs = 5, k = 2, r = 0, level = 7)[["data"]]) -
                 {set.seed(12); unname(generate_artificial_vec(nobs = 5, k = 2, r = 0)[["data"]])},
               matrix(7, 5, 2), ignore_attr = TRUE)

  expect_error(generate_artificial_vec(level = 100), "requires a constant term")
  expect_error(generate_artificial_vec(level = c(1, 2)), "1 or 'k' finite elements")
  expect_error(generate_artificial_vec(const = "restricted", level = NA), "1 or 'k' finite elements")
})
