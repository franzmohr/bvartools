# Forecasts that start from values given in argument 'y_start' of
# add_forecast_input(), which may lie after the end of the estimation sample.
#
# What is checked: that start values equal to the end of the sample reproduce
# the default exactly; that a forecast from start values after the sample
# applies the coefficients to them and to the exogenous values of the new
# periods (the estimand, by hand for an AR(1)); that the deterministic terms,
# seasonal dummies and dummy variables included, are those of the new calendar
# periods; that predict() and add_forecast_errors() date the forecasts by the
# new periods; and what is refused, with which message. The forecast recursion
# itself is tested in test-forecasts.R.

# Quarterly Austrian inflation in percent, an AR model's data
start_inflation <- function(end = c(1998, 1)) {
  stats::window(at_domestic()[, "Dp", drop = FALSE] * 100, end = end)
}

# The change in the oil price in percent, an exogenous variable
start_oil <- function() {
  foreign <- bvartools::at_macrodata[["foreign"]]
  stats::ts(cbind(dpoil = diff(foreign[, "poil"]) * 100), end = stats::end(foreign),
            frequency = 4)
}

# The deterministic columns of forecast regressors
det_columns <- function(object, x) {
  names <- object[["model"]][["deterministic"]]
  det <- x[, ncol(x) - length(names) + seq_along(names), drop = FALSE]
  colnames(det) <- names
  det
}

test_that("start values that end with the sample reproduce the default exactly", {
  model <- fx_var_fitted()
  y <- model[["data"]][["train"]][["y"]]

  default <- add_forecast_input(model, n_ahead = 4)
  started <- add_forecast_input(model, n_ahead = 4, y_start = y)
  expect_identical(started[["data"]][["forecast"]][["x"]], default[["data"]][["forecast"]][["x"]])
  expect_identical(started[["model"]][["h"]], default[["model"]][["h"]])
  expect_equal(started[["model"]][["forecast_start"]], stats::tsp(y)[2] + 1 / 4)

  set.seed(2024)
  default <- add_posterior_forecasts(default)
  set.seed(2024)
  started <- add_posterior_forecasts(started)
  expect_identical(unclass(started[["posterior"]][["forecast"]][["forecasts"]]),
                   unclass(default[["posterior"]][["forecast"]][["forecasts"]]))
  expect_identical(predict(started), predict(default))

  # Also with two lags, an exogenous variable and its lag
  exo <- start_oil()
  y <- var_data()
  model <- create_bvarmodel(y, p = 2, exogen = stats::window(exo, end = stats::end(y)), s = 1,
                            deterministic = "both", iterations = 10, burnin = 5)
  expect_identical(prepare_forecast_input(model, n_ahead = 3, exogen = exo,
                                          y_start = y)[["x"]],
                   prepare_forecast_input(model, n_ahead = 3, exogen = exo)[["x"]])
})

test_that("a forecast from start values after the sample applies the coefficients to them", {
  y <- start_inflation()
  exo <- start_oil()
  model <- create_bvarmodel(y, p = 1, exogen = stats::window(exo, end = stats::end(y)), s = 0,
                            deterministic = "const", iterations = 200, burnin = 50)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  model <- add_posterior_coefficients(add_seed(add_initial_values(model), 31))

  # A nowcast for 1998Q2 far from the last observation, so that a forecast
  # from the end of the sample instead could not pass for one from it
  nowcast <- stats::window(at_domestic()[, "Dp", drop = FALSE] * 100, end = c(1998, 2))
  nowcast[length(nowcast)] <- nowcast[length(nowcast)] + 50
  model <- add_forecast_input(model, n_ahead = 2, exogen = exo, y_start = nowcast)

  # The regressors of 1998Q3: the nowcast, the oil price of 1998Q3, the constant
  x1 <- c(nowcast[length(nowcast)], stats::window(exo, start = c(1998, 3), end = c(1998, 3)), 1)
  expect_equal(model[["data"]][["forecast"]][["x"]][1, ], as.numeric(x1))

  set.seed(32)
  model <- add_posterior_forecasts(model)
  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  fcst <- as.matrix(model[["posterior"]][["forecast"]][["forecasts"]])[, 1]
  precision <- as.numeric(as.matrix(model[["posterior"]][["u_sigma_inv"]][["coeffs"]]))

  # Each draw is its coefficients applied to x1 plus an error drawn with its
  # precision, so the standardised errors are standard normal draws
  z <- (fcst - as.numeric(a %*% x1)) * sqrt(precision)
  expect_lt(abs(mean(z)), 4 / sqrt(length(z)))
  expect_lt(abs(stats::sd(z) - 1), 0.2)
  # The posterior mean forecast is the posterior mean coefficients applied to
  # x1, up to the mean of those errors
  expect_lt(abs(mean(fcst) - sum(colMeans(a) * x1)), 4 * stats::sd(fcst - a %*% x1) / sqrt(length(z)))

  # predict() dates the forecasts by the periods after the nowcast, and shows
  # the nowcast as the history they continue
  pred <- predict(model)
  expect_equal(stats::tsp(pred[["fcst"]])[1:2], c(1998.5, 1998.75))
  expect_equal(stats::tsp(pred[["y"]])[2], 1998.25)
  expect_equal(as.numeric(pred[["y"]][nrow(pred[["y"]]), 1]), as.numeric(nowcast[length(nowcast)]))
  expect_plots(plot(pred))

  # and add_forecast_errors() scores them against the same periods
  realised <- at_domestic()[, "Dp", drop = FALSE] * 100
  scored <- add_forecast_errors(model, test_sample = realised)
  expect_equal(as.numeric(scored[["data"]][["test"]][["y"]]),
               as.numeric(stats::window(realised, start = c(1998, 3), end = c(1998, 4))))
})

test_that("deterministic terms are those of the new calendar periods", {
  y <- start_inflation()
  model <- create_bvarmodel(y, p = 2, deterministic = "both", seasonal = TRUE,
                            iterations = 10, burnin = 5)
  later <- start_inflation(end = c(1998, 3))

  # Starting two quarters later is the default forecast two quarters on
  shifted <- prepare_forecast_input(model, n_ahead = 5, y_start = later)[["x"]]
  default <- prepare_forecast_input(model, n_ahead = 7)[["x"]]
  expect_equal(det_columns(model, shifted), det_columns(model, default)[2 + 1:5, ])

  # and the seasonal dummies of 1998Q4 are those of 1997Q4 in the sample
  train_x <- model[["data"]][["train"]][["x"]]
  season <- grep("season", colnames(train_x), value = TRUE)
  expect_equal(det_columns(model, shifted)[1, season],
               train_x[abs(stats::time(train_x) - 1997.75) < 1e-6, season])
  expect_equal(prepare_forecast_input(model, n_ahead = 5, y_start = later)[["start"]], 1998.75)

  # Monthly data, a sample that ends in April 2018 and a start value for
  # September: the forecast periods are October to December
  set.seed(4711)
  monthly <- stats::ts(cbind(y = sin(1:100 / 3) + stats::rnorm(100, sd = 0.1)),
                       start = c(2010, 1), frequency = 12)
  model <- create_bvarmodel(monthly, p = 1, deterministic = "const", seasonal = TRUE,
                            iterations = 10, burnin = 5)
  start <- stats::ts(cbind(y = 0.5), start = c(2018, 9), frequency = 12)
  shifted <- prepare_forecast_input(model, n_ahead = 3, y_start = start)[["x"]]
  default <- prepare_forecast_input(model, n_ahead = 8)[["x"]]
  expect_equal(shifted[, -1], default[5 + 1:3, -1])
  expect_equal(shifted[1, 1], 0.5)
})

test_that("dummy variables and given deterministic terms follow the new periods", {
  y <- var_data()
  given <- stats::ts(cbind(event = rep(c(0, 1, 0), length.out = nrow(y) + 6)),
                     start = stats::start(y), frequency = 4)
  model <- create_bvarmodel(y, p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- add_dummy_variables(add_dummy_variables(model, step = c(1993, 1)), data = given)
  later <- stats::window(at_domestic(), end = c(1998, 3))
  later <- stats::ts.intersect(y = diff(later[, "y"]), Dp = later[, "Dp"],
                               r = diff(later[, "r"])) * 100

  x <- prepare_forecast_input(model, n_ahead = 3, y_start = later)[["x"]]
  det <- det_columns(model, x)
  expect_equal(unname(det[, "step.1993Q1"]), rep(1, 3))
  expect_equal(unname(det[, "event"]),
               as.numeric(stats::window(given, start = c(1998, 4), end = c(1999, 2))))
  expect_equal(unname(x[1, 1:3]), as.numeric(later[nrow(later), c("y", "Dp", "r")]))

  # Given deterministic terms are read in the new periods
  plain <- create_bvarmodel(y, p = 1, deterministic = "const", iterations = 10, burnin = 5)
  det_given <- stats::ts(cbind(const = 1:6), start = c(1998, 2), frequency = 4)
  x <- prepare_forecast_input(plain, n_ahead = 3, y_start = later, deterministic = det_given)[["x"]]
  expect_equal(x[, 4], 3:5)
})

test_that("start values, exogenous values and deterministic terms that do not fit are refused", {
  y <- start_inflation()
  exo <- start_oil()
  model <- create_bvarmodel(y, p = 2, exogen = stats::window(exo, end = stats::end(y)), s = 1,
                            deterministic = "const", iterations = 10, burnin = 5)
  later <- start_inflation(end = c(1998, 3))

  # Too short for the lags, missing in them, before the end of the sample
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo,
                                  y_start = stats::window(later, start = c(1998, 3))),
               "must contain the last 2 period\\(s\\) before the first forecast period")
  with_na <- later
  with_na[length(with_na)] <- NA
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo, y_start = with_na),
               "has missing values in its last 2 period")
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo,
                                  y_start = start_inflation(end = c(1997, 3))),
               "before the end of the estimation sample")
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo, y_start = as.numeric(later)),
               "must be a time-series object")
  named <- later
  colnames(named) <- "inflation"
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo, y_start = named),
               "must contain the endogenous variables of the model, 'Dp'")
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo,
                                  y_start = stats::ts(as.numeric(later), end = c(1998, 9),
                                                      frequency = 12)),
               "has monthly data, but the model is estimated on quarterly data")

  # Exogenous values that reach the old forecast periods but not the new ones:
  # the lag of 1998Q4 is 1998Q3
  short_exo <- stats::window(exo, end = c(1998, 4))
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = short_exo, y_start = later),
               "must cover the periods from c\\(1998, 3\\) to c\\(1999, 1\\)")
  expect_silent(add_forecast_input(model, n_ahead = 1, exogen = short_exo, y_start = later))

  # Given deterministic terms that do not reach the new forecast periods
  det_given <- stats::ts(cbind(const = rep(1, 3)), start = c(1998, 2), frequency = 4)
  expect_error(add_forecast_input(model, n_ahead = 2, exogen = exo, y_start = later,
                                  deterministic = det_given),
               "Argument 'deterministic' must cover the forecast periods from c\\(1998, 4\\) to c\\(1999, 1\\)")
})

test_that("what cannot start from the start values is refused with a reason", {
  y <- start_inflation()
  later <- start_inflation(end = c(1998, 3))

  # A model that completes a panel not observed whole starts from its draws
  gappy <- y
  gappy[length(gappy) - 3] <- NA
  model <- create_bvarmodel(gappy, p = 1, deterministic = "const", missing = "estimate",
                            iterations = 10, burnin = 5)
  expect_error(add_forecast_input(model, n_ahead = 2, y_start = later),
               "not available for a model estimated with missing = \"estimate\"")

  # Time varying coefficients would not be stepped up to the start values
  model <- create_bvarmodel(y, p = 1, deterministic = "const", tvp = TRUE,
                            iterations = 10, burnin = 5)
  model <- add_forecast_input(model, n_ahead = 2, y_start = later)
  expect_error(add_posterior_forecasts(model), "forecast_states = \"hold\"")

  # A scenario is drawn from the end of the sample
  model <- add_forecast_input(fx_var_fitted(), n_ahead = 2,
                              y_start = stats::ts(matrix(0, 1, 3, dimnames = list(NULL, c("y", "Dp", "r"))),
                                                  start = c(1998, 3), frequency = 4))
  expect_equal(model[["model"]][["forecast_start"]], 1998.75)
  scenario <- data.frame(period = 1, variable = "y", value = 1)
  expect_error(add_posterior_forecasts(model, scenario = scenario), "a scenario is not available")

  # So is a score; a VEC model and the windows of an expanding window forecast
  # from the end of their own samples
  forecasted <- add_posterior_forecasts(model)
  expect_error(add_predictive_loglik(forecasted, test_sample = var_data()),
               "a predictive log likelihood is not available")
  expect_error(add_forecast_input(fx_vec_model(), n_ahead = 2, y_start = vec_data()),
               "not available for a 'bvecmodel'")
  windows <- use_expanding_window(fx_var_model(), start = 1996)
  expect_error(add_forecast_input(windows, n_ahead = 2, y_start = var_data()),
               "not available for an 'expandingwindow' object")

  # Without start values the model forecasts from the sample again
  expect_null(add_forecast_input(model, n_ahead = 2)[["model"]][["forecast_start"]])
})
