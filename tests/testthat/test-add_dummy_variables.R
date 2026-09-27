# Dummy variables added to the deterministic terms of a model.
#
# What is checked: that a dummy lands in the period it was asked for whatever
# sample the lags leave, that the SUR form around it is rebuilt without losing
# the columns of a structural or a VEC model, that an impulse dummy takes up its
# observation (the estimand), that a forecast continues each dummy by its rule,
# that the windows of an expanding window leave out a dummy they do not reach,
# that the rule survives a file, and what is refused. The priors and samplers
# treat a dummy as any deterministic term and are tested elsewhere.

dummy_model <- function(p = 1, ...) {
  create_bvarmodel(var_data(), p = p, deterministic = "const",
                   iterations = fx_iterations, burnin = fx_burnin, ...)
}

x_at <- function(object, name) {
  x <- object[["data"]][["train"]][["x"]]
  stats::time(x)[x[, name] == 1]
}

sur_form <- function(x, k) {
  z <- kronecker(as.matrix(x), diag(1, k))
  dimnames(z) <- NULL
  z
}

test_that("impulse and step dummies are built in the periods asked for", {
  object <- add_dummy_variables(dummy_model(), impulse = list(c(1990, 1), 1991.5),
                                step = c(1993, 1))

  expect_identical(object[["model"]][["deterministic"]],
                   c("const", "impulse.1990Q1", "impulse.1991Q3", "step.1993Q1"))
  expect_identical(object[["model"]][["n"]], 4L)
  expect_identical(object[["model"]][["dummy_variables"]][["type"]],
                   c("impulse", "impulse", "step"))
  expect_equal(as.numeric(x_at(object, "impulse.1990Q1")), 1990)
  expect_equal(as.numeric(x_at(object, "impulse.1991Q3")), 1991.5)
  expect_equal(min(x_at(object, "step.1993Q1")), 1993)
  expect_equal(max(x_at(object, "step.1993Q1")), max(stats::time(object[["data"]][["train"]][["y"]])))

  # The SUR form is the compact regressors, dummies included.
  expect_identical(object[["data"]][["train"]][["z"]], sur_form(object[["data"]][["train"]][["x"]], 3))
})

test_that("a dummy lands in the same period in every model of a list", {
  models <- add_dummy_variables(dummy_model(p = 1:3), impulse = c(1990, 1))

  expect_s3_class(models, "modellist")
  for (object in models) {
    expect_equal(as.numeric(x_at(object, "impulse.1990Q1")), 1990)
  }
})

test_that("a series given in 'data' takes the place a built-in term has", {
  both <- create_bvarmodel(var_data(), p = 1, deterministic = "both",
                           iterations = fx_iterations, burnin = fx_burnin)
  linear <- both[["data"]][["train"]][["x"]][, "trend", drop = FALSE]
  colnames(linear) <- "linear"

  object <- add_dummy_variables(dummy_model(), data = linear)

  expect_equal(unname(as.matrix(object[["data"]][["train"]][["x"]])),
               unname(as.matrix(both[["data"]][["train"]][["x"]])))
  expect_identical(object[["data"]][["train"]][["z"]], both[["data"]][["train"]][["z"]])
})

test_that("the contemporaneous columns of a structural model are kept", {
  object <- create_bvarmodel(var_data(), p = 1, structural = TRUE, error = "gamma",
                             iterations = fx_iterations, burnin = fx_burnin)
  z_old <- object[["data"]][["train"]][["z"]]
  n_old <- ncol(object[["data"]][["train"]][["x"]])

  object <- add_dummy_variables(object, impulse = c(1990, 1))
  z <- object[["data"]][["train"]][["z"]]

  expect_identical(z[, seq_len(3 * (n_old + 1))], sur_form(object[["data"]][["train"]][["x"]], 3))
  expect_identical(z[, -seq_len(3 * (n_old + 1))], z_old[, -seq_len(3 * n_old)])
})

test_that("an impulse dummy with a flat prior takes up its observation", {
  object <- add_dummy_variables(dummy_model(), impulse = c(1990, 1))
  object[["model"]][["iterations"]] <- 2000L
  object[["model"]][["burnin"]] <- 200L
  object <- add_priors(object, coef = list(v_i = 1 / 9, v_i_det = 0),
                       sigma = list(df = "k", scale = 1))
  object <- add_seed(object, 11)
  object <- add_posterior_coefficients(add_initial_values(object))

  y <- object[["data"]][["train"]][["y"]]
  x <- as.matrix(object[["data"]][["train"]][["x"]])
  a <- matrix(colMeans(object[["posterior"]][["a"]][["coeffs"]]), 3)
  residuals <- as.matrix(y) - x %*% t(a)
  t_star <- which(abs(stats::time(y) - 1990) < 1e-6)

  expect_lt(max(abs(residuals[t_star, ])), 0.1 * min(apply(residuals[-t_star, ], 2, stats::sd)))
})

test_that("a forecast continues impulse dummies with zero and step dummies with one", {
  object <- add_dummy_variables(dummy_model(), impulse = c(1990, 1), step = c(1993, 1))
  object <- add_forecast_input(object, n_ahead = 3)
  x <- object[["data"]][["forecast"]][["x"]]
  names <- object[["model"]][["deterministic"]]
  det <- x[, ncol(x) - length(names) + seq_along(names), drop = FALSE]
  colnames(det) <- names

  expect_equal(unname(det[, "const"]), rep(1, 3))
  expect_equal(unname(det[, "impulse.1990Q1"]), rep(0, 3))
  expect_equal(unname(det[, "step.1993Q1"]), rep(1, 3))
})

test_that("a series in 'data' is continued with the values it was given", {
  y <- var_data()
  given <- stats::ts(cbind(event = rep(c(0, 1), length.out = nrow(y) + 3)),
                     start = stats::start(y), frequency = 4)
  object <- add_dummy_variables(dummy_model(), data = given)
  object <- add_forecast_input(object, n_ahead = 3)
  x <- object[["data"]][["forecast"]][["x"]]
  expected <- as.numeric(stats::window(given, start = c(1998, 2)))
  expect_equal(x[, ncol(x)], expected)

  short <- stats::window(given, end = c(1998, 1))
  object <- add_dummy_variables(dummy_model(), data = short)
  expect_error(add_forecast_input(object, n_ahead = 3), "without values for the forecast periods")
})

test_that("what cannot be estimated is refused with a reason", {
  object <- dummy_model()

  expect_error(add_dummy_variables(object), "at least one dummy variable")
  expect_error(add_dummy_variables(object, impulse = c(2010, 1)), "zero in every period")
  expect_error(add_dummy_variables(object, step = c(1970, 1)), "is the constant again")
  expect_error(add_dummy_variables(object, impulse = c(1990, 5)), "does not exist")
  expect_error(add_dummy_variables(object, impulse = list(c(1990, 1), 1990)), "specified twice")
  expect_error(add_dummy_variables(add_dummy_variables(object, impulse = c(1990, 1)),
                                   impulse = c(1990, 1)),
               "already has")

  named <- stats::ts(cbind(trend = seq_len(nrow(var_data()))), start = stats::start(var_data()),
                     frequency = 4)
  expect_error(add_dummy_variables(object, data = named), "are taken by the deterministic terms")
  expect_error(add_dummy_variables(object, data = stats::ts(1:100, start = 1979, frequency = 4)),
               "must be named")

  with_priors <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                            sigma = list(df = "k", scale = 1))
  expect_error(add_dummy_variables(with_priors, impulse = c(1990, 1)), "before add_priors")
})

test_that("windows that end before a dummy's period leave it out", {
  start <- c(1989, 3)
  by_windows <- add_dummy_variables(use_expanding_window(dummy_model(), start = start),
                                    impulse = c(1990, 1))
  by_model <- use_expanding_window(add_dummy_variables(dummy_model(), impulse = c(1990, 1)),
                                   start = start)

  for (windows in list(by_windows, by_model)) {
    ends <- vapply(windows, function(w) max(stats::time(w[["data"]][["train"]][["y"]])), 0)
    has <- vapply(windows, function(w) "impulse.1990Q1" %in% w[["model"]][["deterministic"]], TRUE)
    expect_identical(has, ends >= 1990)
    for (w in windows) {
      expect_identical(ncol(w[["data"]][["train"]][["x"]]), 3L + w[["model"]][["n"]])
      expect_identical(w[["data"]][["train"]][["z"]], sur_form(w[["data"]][["train"]][["x"]], 3))
    }
  }
  expect_equal(by_windows, by_model)

  with_priors <- add_priors(add_dummy_variables(dummy_model(), impulse = c(1990, 1)),
                            coef = list(v_i = 1, v_i_det = 0.1),
                            sigma = list(df = "k", scale = 1))
  expect_error(suppressWarnings(use_expanding_window(with_priors, start = start)),
               "before add_priors")
})

test_that("a VEC model keeps its cointegration columns and passes the dummies on", {
  object <- create_bvecmodel(at_data(), p = 2, r = 1, const = "unrestricted",
                             iterations = fx_iterations, burnin = fx_burnin)
  z_old <- object[["data"]][["train"]][["z"]]
  object <- add_dummy_variables(object, impulse = c(2008, 4))
  z <- object[["data"]][["train"]][["z"]]

  expect_identical(ncol(z), ncol(z_old) + 3L)
  expect_identical(z[, 1:3], z_old[, 1:3])
  expect_identical(z[, -(1:3)], sur_form(object[["data"]][["train"]][["x"]], 3))

  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = "k", scale = 1))
  object <- add_posterior_coefficients(add_initial_values(object))
  var <- vec_to_var(object)
  expect_identical(var[["model"]][["dummy_variables"]], object[["model"]][["dummy_variables"]])
  expect_true("impulse.2008Q4" %in% var[["model"]][["deterministic"]])

  object <- add_posterior_forecasts(add_forecast_input(object, n_ahead = 2))
  expect_false(anyNA(object[["posterior"]][["forecast"]][["forecasts"]]))
})

test_that("a VEC model continues a series in 'data' through its VAR representation", {
  object <- create_bvecmodel(at_data(), p = 2, r = 1, const = "unrestricted",
                             iterations = fx_iterations, burnin = fx_burnin)
  y <- object[["data"]][["train"]][["y"]]
  given <- stats::ts(cbind(event = rep(c(0, 1, 0), length.out = nrow(y) + 2)),
                     start = stats::start(y), frequency = 4)
  object <- add_dummy_variables(object, data = given)
  object <- add_forecast_input(object, n_ahead = 2)
  x <- object[["data"]][["forecast"]][["x"]]

  expect_equal(x[, ncol(x)], as.numeric(given)[nrow(y) + 1:2])
})

test_that("the rule a dummy is continued by survives a file", {
  skip_if_not_installed("hdf5r")
  object <- add_dummy_variables(dummy_model(), impulse = c(1990, 1), step = c(1993, 1))
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path))

  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(restored[["model"]][["dummy_variables"]], object[["model"]][["dummy_variables"]])
  expect_equal(add_forecast_input(restored, n_ahead = 2)[["data"]][["forecast"]],
               add_forecast_input(object, n_ahead = 2)[["data"]][["forecast"]])
})
