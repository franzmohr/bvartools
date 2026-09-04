test_that("add_forecast_input prepares the regressors of the forecast periods", {
  model <- add_forecast_input(fx_var_fitted(), n_ahead = 5)
  spec <- model[["model"]]

  expect_equal(spec[["h"]], 5)
  expect_false(is.null(model[["data"]][["forecast"]][["z"]]))
  # The forecast regressors are in SUR form: one row block per period.
  expect_identical(nrow(model[["data"]][["forecast"]][["z"]]),
                   as.integer(5 * spec[["k"]]))
  expect_identical(ncol(model[["data"]][["forecast"]][["z"]]),
                   ncol(model[["data"]][["train"]][["z"]]))
})

test_that("prepare_forecast_input returns the horizon and its regressors", {
  input <- prepare_forecast_input(fx_var_fitted(), n_ahead = 3)

  expect_named(input, c("h", "z"))
  expect_identical(input[["h"]], 3L)
  expect_identical(nrow(input[["z"]]),
                   as.integer(3 * fx_var_fitted()[["model"]][["k"]]))
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
  train <- stats::window(data, end = c(1977, 4))
  test <- stats::window(data, start = c(1978, 1))

  model <- create_bvarmodel(train, p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  model <- add_initial_values(model)
  set.seed(31)
  model <- add_posterior_coefficients(model)
  model <- add_forecast_input(model, n_ahead = 4)
  model <- add_posterior_forecasts(model)
  model <- add_forecast_errors(model, test_sample = test)

  errors <- get_forecast_errors(model)

  expect_false(is.null(errors))
  expect_true(length(errors) > 0)
})
