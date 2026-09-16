test_that("add_forecast_input prepares the regressors of the forecast periods", {
  model <- add_forecast_input(fx_var_fitted(), n_ahead = 5)
  spec <- model[["model"]]

  expect_equal(spec[["h"]], 5)
  expect_false(is.null(model[["data"]][["forecast"]][["x"]]))
  # The forecast regressors are compact: one row per period, one column per
  # regressor, so k times narrower and k times shorter than the training SUR
  # matrix beside them.
  expect_identical(nrow(model[["data"]][["forecast"]][["x"]]), 5L)
  expect_identical(ncol(model[["data"]][["forecast"]][["x"]]) * spec[["k"]],
                   ncol(model[["data"]][["train"]][["z"]]))
})

test_that("prepare_forecast_input returns the horizon and its regressors", {
  input <- prepare_forecast_input(fx_var_fitted(), n_ahead = 3)

  expect_named(input, c("h", "x"))
  expect_identical(input[["h"]], 3L)
  expect_identical(nrow(input[["x"]]), 3L)
})

test_that("a forecast input in the old SUR layout is still accepted", {
  # What an object fitted before the compact layout carries. The C++ side
  # compacts it on the way in, so the draws have to match what the compact
  # spelling of the same regressors produces.
  model <- add_forecast_input(fx_var_fitted(), n_ahead = 4)
  k <- model[["model"]][["k"]]

  legacy <- model
  legacy[["data"]][["forecast"]] <- list("z" = kronecker(model[["data"]][["forecast"]][["x"]],
                                                        diag(1, k)))

  set.seed(7357)
  compact <- add_posterior_forecasts(model)
  set.seed(7357)
  old <- add_posterior_forecasts(legacy)

  expect_equal(unname(as.matrix(old[["posterior"]][["forecast"]])),
               unname(as.matrix(compact[["posterior"]][["forecast"]])))
})

test_that("forecast draws are stored for every period and variable", {
  model <- fx_var_forecast()
  spec <- model[["model"]]
  draws <- model[["posterior"]][["forecast"]]

  expect_s3_class(draws, "mcmc")
  expect_identical(dim(draws), c(fx_iterations, as.integer(5 * spec[["k"]])))
  expect_true(all(is.finite(draws)))
})

test_that("predict reshapes the draws into horizon, variable and iteration", {
  model <- fx_var_forecast()
  forecast <- stats::predict(model, n_ahead = 5)

  expect_s3_class(forecast, "bvarprd")
  expect_named(forecast, c("fcst", "y"))
  expect_identical(dim(forecast[["fcst"]]),
                   c(5L, model[["model"]][["k"]], fx_iterations))
  expect_identical(dimnames(forecast[["fcst"]])[[2]],
                   model[["model"]][["endogen"]])

  # Each draw is the corresponding row of the stored draws, read variable by
  # variable within a period.
  expected <- t(matrix(model[["posterior"]][["forecast"]][1, ],
                       model[["model"]][["k"]]))
  expect_equal(unname(forecast[["fcst"]][, , 1]), unname(expected))
})

test_that("the forecast continues the sample in time", {
  model <- fx_var_forecast()
  forecast <- stats::predict(model, n_ahead = 5)

  y <- model[["data"]][["train"]][["y"]]
  frequency <- stats::frequency(y)
  last_period <- stats::tsp(y)[2]

  expect_equal(stats::tsp(forecast[["fcst"]])[1], last_period + 1 / frequency)
  expect_equal(stats::tsp(forecast[["fcst"]])[3], frequency)
})

test_that("predict warns when more periods are requested than simulated", {
  expect_warning(stats::predict(fx_var_forecast(), n_ahead = 20),
                 "larger than the value")
})

test_that("a shorter horizon can be requested", {
  forecast <- stats::predict(fx_var_forecast(), n_ahead = 2)
  expect_identical(dim(forecast[["fcst"]])[1], 2L)
})

test_that("predict needs the forecast input and the simulated draws", {
  expect_error(stats::predict(fx_var_fitted()), "Missing specification of h")
  expect_error(
    stats::predict(add_forecast_input(fx_var_fitted(), n_ahead = 3)),
    "Missing element"
  )
})

test_that("forecast errors are computed against a test sample", {
  data <- var_data()
  train <- stats::window(data, end = c(1997, 1))
  test <- stats::window(data, start = c(1997, 2))

  model <- create_bvarmodel(train, p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  set.seed(31)
  model <- add_initial_values(model)
  model <- add_posterior_coefficients(model)
  model <- add_forecast_input(model, n_ahead = 4)
  model <- add_posterior_forecasts(model)
  model <- add_forecast_errors(model, test_sample = test)

  errors <- get_forecast_errors(model)

  expect_false(is.null(errors))
  expect_true(length(errors) > 0)
})

test_that("a time varying model simulates its states forward unless told to hold them", {
  for (error in c("sv", "gamma")) {
    fitted <- add_forecast_input(fx_var_tvp_fitted(error), n_ahead = 6)
    k <- fitted[["model"]][["k"]]

    if (error == "sv") {
      # The step the volatility is simulated forward by, one per variable and draw.
      sigma <- fitted[["posterior"]][["u_sigma_inv"]][["sigma"]]
      expect_s3_class(sigma, "mcmc")
      expect_identical(dim(sigma), c(nrow(fitted[["posterior"]][["a"]][["coeffs"]]), k))
    }

    set.seed(2718)
    simulated <- add_posterior_forecasts(fitted)
    set.seed(2718)
    held <- add_posterior_forecasts(fitted, forecast_states = "hold")

    expect_true(all(is.finite(simulated[["posterior"]][["forecast"]])))
    expect_identical(held[["model"]][["forecast_states"]], "hold")
    # From the same seed, the drift is the one thing separating the two.
    expect_false(isTRUE(all.equal(unclass(simulated[["posterior"]][["forecast"]]),
                                  unclass(held[["posterior"]][["forecast"]]))))
  }

  expect_error(add_posterior_forecasts(add_forecast_input(fx_var_tvp_fitted("gamma"), n_ahead = 2),
                                       forecast_states = "sideways"))
})

test_that("a stochastic volatility VEC model keeps the variance of its log-volatility innovations", {
  fitted <- fx_vec_tvp_fitted("sv")
  sigma <- fitted[["posterior"]][["u_sigma_inv"]][["sigma"]]

  # One per variable and draw, a chain like the other blocks, as the VAR models
  # keep it -- even though a VEC model's forecast still holds its volatility.
  expect_s3_class(sigma, "mcmc")
  expect_identical(dim(sigma), c(nrow(fitted[["posterior"]][["a"]][["coeffs"]]),
                                 fitted[["model"]][["k"]]))
  expect_true(all(sigma > 0))
})

test_that("a stochastic volatility posterior without volatility steps forecasts only when held", {
  old <- add_forecast_input(fx_var_tvp_fitted("sv"), n_ahead = 3)
  old[["posterior"]][["u_sigma_inv"]][["sigma"]] <- NULL

  expect_error(add_posterior_forecasts(old), "innovation variances")
  expect_no_error(add_posterior_forecasts(old, forecast_states = "hold"))
})

test_that("predict returns every simulated period by default", {
  expect_no_warning(forecast <- stats::predict(fx_var_forecast()))
  expect_identical(dim(forecast[["fcst"]])[1], fx_var_forecast()[["model"]][["h"]])
})
