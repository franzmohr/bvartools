test_that("vec_to_var produces a VAR model of the same size", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)

  expect_s3_class(var, "bvarmodel")
  expect_identical(var[["model"]][["type"]], "VAR")
  expect_identical(var[["model"]][["k"]], vec[["model"]][["k"]])
  expect_identical(var[["model"]][["endogen"]], vec[["model"]][["endogen"]])
  # p is the lag order of the model in levels in both parameterisations.
  expect_equal(var[["model"]][["p"]], vec[["model"]][["p"]])
})

test_that("the converted model carries the posterior draws over", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)
  spec <- var[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_identical(nrow(var[["posterior"]][["a"]][["coeffs"]]), fx_iterations)
  expect_identical(ncol(var[["posterior"]][["a"]][["coeffs"]]),
                   as.integer(n_coeffs))
  # The error precision is unchanged by the reparameterisation.
  expect_equal(var[["posterior"]][["u_sigma_inv"]][["coeffs"]],
               vec[["posterior"]][["u_sigma_inv"]][["coeffs"]])
})

test_that("the level coefficients follow from Pi and the short-run terms", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)
  k <- vec[["model"]][["k"]]
  r <- vec[["model"]][["rank"]]

  draw <- 1L
  alpha <- matrix(vec[["posterior"]][["a"]][["coeffs"]][draw, seq_len(k * r)], k)
  beta <- matrix(vec[["posterior"]][["beta"]][["coeffs"]][draw, ],
                 vec[["model"]][["k_beta"]])

  # With p = 1 there are no lagged differences, so A_1 = Pi + I.
  expected <- alpha %*% t(beta) + diag(1, k)
  actual <- matrix(var[["posterior"]][["a"]][["coeffs"]][draw, seq_len(k^2)], k)
  expect_equal(actual, expected)
})

test_that("the converted model is usable by the application functions", {
  var <- vec_to_var(fx_vec_fitted())

  # The deterministic terms of the VEC keep their names across the conversion.
  expect_identical(var[["model"]][["deterministic"]], "const")

  expect_s3_class(irf(var, impulse = "lr", response = "Dp", n_ahead = 3),
                  "bvarirf")
  expect_s3_class(fevd(var, response = "Dp", n_ahead = 3), "bvarfevd")
  expect_no_error(summary(var))
})

test_that("a modellist of VEC models is converted element by element", {
  models <- create_bvecmodel(vec_data(), p = 1, r = 1:2,
                             const = "unrestricted",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 1, v_i_det = 1 / 10),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = "k", scale = 1))
  set.seed(17)
  models <- add_initial_values(models)
  models <- add_posterior_coefficients(models)
  converted <- vec_to_var(models)

  expect_s3_class(converted, "modellist")
  expect_length(converted, 2)
  expect_true(all(vapply(converted, inherits, logical(1), "bvarmodel")))
})

test_that("a time varying VEC model is converted period by period", {
  vec <- fx_vec_tvp_fitted("sv")
  var <- vec_to_var(vec)
  k <- vec[["model"]][["k"]]
  r <- vec[["model"]][["rank"]]
  tt <- nrow(vec[["data"]][["train"]][["y"]])
  n_vec <- ncol(vec[["data"]][["train"]][["z"]])
  n_var <- ncol(var[["data"]][["train"]][["z"]])
  k_beta <- ncol(vec[["data"]][["train"]][["w"]])

  expect_s3_class(var, "bvarmodel")
  expect_true(var[["model"]][["tvp"]])
  expect_identical(var[["model"]][["algorithm"]], "VarTvpStochvol")
  expect_identical(dim(var[["posterior"]][["a"]][["coeffs"]]),
                   c(nrow(vec[["posterior"]][["a"]][["coeffs"]]), as.integer(n_var * tt)))
  # The drift of the VEC coefficients and the cointegration space have no
  # counterpart in levels; the error term is shared.
  expect_null(var[["posterior"]][["beta"]])
  expect_null(var[["posterior"]][["a"]][["sigma"]])
  expect_equal(var[["posterior"]][["u_sigma_inv"]][["coeffs"]],
               vec[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  # With p = 2, A_1 = I + Pi_t + Gamma_1 and A_2 = -Gamma_1 in every period.
  draw <- 1L
  for (period in c(1L, tt)) {
    a_vec <- vec[["posterior"]][["a"]][["coeffs"]][draw, (period - 1) * n_vec + seq_len(n_vec)]
    alpha <- matrix(a_vec[seq_len(k * r)], k)
    gamma <- matrix(a_vec[k * r + seq_len(k^2)], k)
    beta <- matrix(vec[["posterior"]][["beta"]][["coeffs"]][draw, (period - 1) * k_beta * r +
                                                              seq_len(k_beta * r)], k_beta)
    a_var <- var[["posterior"]][["a"]][["coeffs"]][draw, (period - 1) * n_var + seq_len(n_var)]

    expect_equal(matrix(a_var[seq_len(k^2)], k),
                 diag(1, k) + alpha %*% t(beta[seq_len(k), , drop = FALSE]) + gamma)
    expect_equal(matrix(a_var[k^2 + seq_len(k^2)], k), -gamma)
  }
})

test_that("a converted time varying VEC model forecasts and responds", {
  var <- vec_to_var(fx_vec_tvp_fitted("gamma"))

  expect_no_error(summary(var, period = 1))
  expect_s3_class(irf(var, impulse = "lr", response = "Dp", n_ahead = 3), "bvarirf")
  expect_s3_class(fevd(var, response = "Dp", n_ahead = 3, period = 1), "bvarfevd")
  expect_plots(plot(var))

  forecast <- add_posterior_forecasts(add_forecast_input(var, n_ahead = 2))
  expect_identical(nrow(forecast[["posterior"]][["forecast"]]), as.integer(fx_iterations))
  expect_identical(ncol(forecast[["posterior"]][["forecast"]]),
                   as.integer(2 * var[["model"]][["k"]]))
})

test_that("a VEC model with stochastic volatility forecasts from its last period", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                            error = "sv", iterations = fx_iterations, burnin = fx_burnin)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = tvp_sigma_prior("sv"))
  set.seed(204)
  model <- add_initial_values(model)
  vec <- add_posterior_coefficients(model)

  # The coefficients are constant, the error precision is a path. The
  # transformation used to be handed the whole path and refused it.
  var <- vec_to_var(vec)
  k <- vec[["model"]][["k"]]
  r <- vec[["model"]][["rank"]]
  kk <- k * k
  tt <- nrow(vec[["data"]][["train"]][["y"]])

  expect_identical(var[["model"]][["algorithm"]], "VarNormalStochvol")
  expect_equal(var[["posterior"]][["u_sigma_inv"]][["coeffs"]],
               vec[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  expect_identical(ncol(var[["posterior"]][["u_sigma_inv"]][["coeffs"]]), as.integer(kk * tt))

  # With p = 2, A_1 = I + Pi + Gamma_1 and A_2 = -Gamma_1.
  draw <- 1L
  a_vec <- vec[["posterior"]][["a"]][["coeffs"]][draw, ]
  alpha <- matrix(a_vec[seq_len(k * r)], k)
  gamma <- matrix(a_vec[k * r + seq_len(kk)], k)
  beta <- matrix(vec[["posterior"]][["beta"]][["coeffs"]][draw, ], vec[["model"]][["k_beta"]])
  a_var <- var[["posterior"]][["a"]][["coeffs"]][draw, ]
  expect_equal(matrix(a_var[seq_len(kk)], k),
               diag(1, k) + alpha %*% t(beta[seq_len(k), , drop = FALSE]) + gamma)
  expect_equal(matrix(a_var[kk + seq_len(kk)], k), -gamma)

  var <- add_forecast_input(var, n_ahead = 2)

  # The coefficients are constant and the volatility is the VEC model's own
  # random walk, whose innovation variances the posterior carries: the VAR
  # representation simulates it forward rather than holding it.
  expect_null(var[["model"]][["forecast_states"]])
  expect_false(is.null(var[["posterior"]][["u_sigma_inv"]][["sigma"]]))

  forecast_with <- function(object, states) {
    set.seed(1)
    unclass(add_posterior_forecasts(object, forecast_states = states)[["posterior"]][["forecast"]])
  }
  simulated <- forecast_with(var, "simulate")
  expect_identical(dim(simulated), c(as.integer(fx_iterations), as.integer(2 * k)))
  expect_true(all(is.finite(simulated)))
  expect_false(isTRUE(all.equal(simulated, forecast_with(var, "hold"),
                                check.attributes = FALSE)))

  # A held forecast draws its errors under the precision of the last in-sample
  # period: a path that holds that period throughout gives the same forecast,
  # one that holds the first period does not. The variants are built before
  # any forecast is, so that each starts from a posterior without one.
  hold_period <- function(object, period) {
    path <- as.matrix(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    held <- path[, rep((period - 1) * kk + seq_len(kk), tt), drop = FALSE]
    object[["posterior"]][["u_sigma_inv"]][["coeffs"]] <- coda::mcmc(held)
    object
  }
  held_last <- hold_period(var, tt)
  held_first <- hold_period(var, 1L)

  forecast <- forecast_with(var, "hold")
  expect_equal(forecast_with(held_last, "hold"), forecast, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(forecast_with(held_first, "hold"), forecast,
                                check.attributes = FALSE)))

  # A posterior without the innovation variances -- one estimated before they
  # were kept -- converts to a VAR representation that holds the volatility.
  old <- vec
  old[["posterior"]][["u_sigma_inv"]][["sigma"]] <- NULL
  expect_identical(vec_to_var(old)[["model"]][["forecast_states"]], "hold")
})

test_that("a converted time varying VEC model holds its states", {
  for (error in c("gamma", "sv")) {
    expect_identical(vec_to_var(fx_vec_tvp_fitted(error))[["model"]][["forecast_states"]], "hold")
  }
})

test_that("an expanding window of VEC models can be evaluated out of sample", {
  full <- vec_data()
  train <- stats::window(full, end = c(1994, 4))
  model <- create_bvecmodel(train, p = 2, r = 1, const = "unrestricted",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1))
  windows <- use_expanding_window(model, start = c(1994, 2))
  set.seed(23)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)

  converted <- vec_to_var(windows)
  expect_s3_class(converted, "expandingwindow")
  expect_true(all(vapply(converted, inherits, logical(1), "bvarmodel")))

  converted <- add_forecast_input(converted, n_ahead = 2)
  converted <- add_posterior_forecasts(converted)
  converted <- add_forecast_errors(converted, test_sample = full)
  criteria <- selection_criteria(converted)
  expect_true(all(c("FE", "AFE", "RSFE") %in% names(criteria)))
})

test_that("impulse responses and decompositions point a VEC model to its VAR representation", {
  # Forecasts are not among them: a 'bvecmodel' is forecast directly, see
  # test-forecast_bvecmodel.R.
  model <- fx_at_vec_tvp()
  calls <- list("irf" = function() irf(model),
                "fevd" = function() fevd(model),
                "spillover" = function() spillover(model))

  for (name in names(calls)) {
    expect_error(calls[[name]](), paste0("'", name, "' does not work directly on a 'bvecmodel'"),
                 label = name)
  }
})
