# The reference model of the comparisons below. The training samples of its
# windows end in 1976Q4 to 1977Q4, the test data cover 1978.
ext_reference <- function() {
  fx_expanding_forecast()
}

ext_ends <- function() {
  unname(sapply(ext_reference(), function(x) {
    stats::tsp(x[["data"]][["train"]][["y"]])[2]
  }))
}

# Point forecasts of the periods that follow a training sample. A publication is
# dated one quarter after the end of the training sample it belongs to, plus a
# few days, which is the situation that argument 'data_lag' describes.
ext_forecasts <- function(value = 0, offset = 0.05) {
  full <- stats::window(diff(log(bvartools::e1)) * 100, end = c(1978, 4))
  result <- NULL
  for (end in ext_ends()) {
    for (h in 1:2) {
      for (variable in dimnames(full)[[2]]) {
        result <- rbind(result,
                        data.frame(origin = end + 0.25 + offset,
                                   period = end + h / 4,
                                   variable = variable,
                                   value = ifelse(is.function(value),
                                                  value(end + h / 4, variable),
                                                  value),
                                   stringsAsFactors = FALSE))
      }
    }
  }
  result
}

test_that("each publication becomes one window of the training sample before it", {
  external <- create_external_forecast(ext_forecasts(), ext_reference(),
                                       n_ahead = 2, data_lag = 1)

  expect_s3_class(external, "externalforecast")
  expect_s3_class(external, "expandingwindow")
  expect_length(external, length(ext_reference()))

  # A forecast published one quarter after the end of a training sample is
  # matched to that training sample, not to a later one
  ends <- sapply(external, function(x) {stats::tsp(x[["data"]][["train"]][["y"]])[2]})
  expect_equal(unname(ends), ext_ends())
})

test_that("the periods of a publication are translated into forecast horizons", {
  external <- create_external_forecast(ext_forecasts(value = 1), ext_reference(),
                                       n_ahead = 2, data_lag = 1)
  forecast <- external[[1]][["posterior"]][["forecast"]]

  # One row, because external forecasts are point forecasts, and one column per
  # endogenous variable within each of the two forecast horizons
  expect_identical(dim(forecast), c(1L, 6L))
  expect_true(all(forecast == 1))
})

test_that("perfect foresight leaves no forecast errors", {
  full <- stats::window(diff(log(bvartools::e1)) * 100, end = c(1978, 4))
  realised <- function(period, variable) {
    as.numeric(stats::window(full[, variable], start = period, end = period))
  }

  external <- create_external_forecast(ext_forecasts(value = realised),
                                       ext_reference(), n_ahead = 2, data_lag = 1)
  external <- add_forecast_errors(external, test_sample = full)

  errors <- unlist(lapply(external, get_forecast_errors))
  expect_true(length(errors) > 0)
  expect_equal(errors, rep(0, length(errors)), ignore_attr = TRUE)
})

test_that("a forecast horizon beyond 'n_ahead' is dropped", {
  forecasts <- ext_forecasts(value = 1)
  # A third period, which is not covered by the forecast horizon of the models
  extra <- forecasts[forecasts[, "period"] == forecasts[1, "origin"] - 0.05 + 0.25, ]
  extra[, "period"] <- extra[, "period"] + 0.75
  external <- create_external_forecast(rbind(forecasts, extra), ext_reference(),
                                       n_ahead = 2, data_lag = 1)

  expect_identical(ncol(external[[1]][["posterior"]][["forecast"]]), 6L)
})

test_that("variables that are not part of the models are reported and dropped", {
  forecasts <- ext_forecasts(value = 1)
  other <- forecasts[1:2, ]
  other[, "variable"] <- "gdp"

  expect_warning(external <- create_external_forecast(rbind(forecasts, other),
                                                      ext_reference(),
                                                      n_ahead = 2, data_lag = 1),
                 "gdp")
  expect_identical(external[[1]][["model"]][["endogen"]],
                   dimnames(ext_reference()[[1]][["data"]][["train"]][["y"]])[[2]])
})

test_that("only one publication per training sample is used", {
  early <- ext_forecasts(value = 1, offset = 0.02)
  late <- ext_forecasts(value = 2, offset = 0.2)
  forecasts <- rbind(early, late)

  last <- create_external_forecast(forecasts, ext_reference(), n_ahead = 2,
                                   data_lag = 1, select = "last")
  first <- create_external_forecast(forecasts, ext_reference(), n_ahead = 2,
                                    data_lag = 1, select = "first")

  expect_length(last, length(ext_ends()))
  expect_length(first, length(ext_ends()))
  expect_true(all(last[[1]][["posterior"]][["forecast"]] == 2))
  expect_true(all(first[[1]][["posterior"]][["forecast"]] == 1))
})

test_that("one object per forecaster is produced", {
  forecasts <- ext_forecasts(value = 1)
  forecasts[["forecaster"]] <- "A"
  other <- forecasts
  other[["forecaster"]] <- "B"

  external <- create_external_forecast(rbind(forecasts, other), ext_reference(),
                                       n_ahead = 2, by = "forecaster", data_lag = 1)

  expect_s3_class(external, "modellist")
  expect_named(external, c("A", "B"))
  expect_s3_class(external[[1]], "externalforecast")
})

test_that("annual forecasts of the publication year are one-step ahead forecasts", {
  annual <- stats::ts(matrix(1:20, 20, 1, dimnames = list(NULL, "gdp")),
                      start = 2001, frequency = 1)
  model <- create_bvarmodel(annual, p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- use_expanding_window(model, start = 2018)

  # Publication dates, as they are provided by external sources
  forecasts <- data.frame(origin = c("2018-11-15", "2018-11-15", "2019-11-15"),
                          period = c(2018, 2019, 2019),
                          variable = "gdp",
                          value = c(1, 2, 3),
                          stringsAsFactors = FALSE)

  external <- create_external_forecast(forecasts, model, n_ahead = 2, data_lag = 1)

  # The data of 2018 were not available in November 2018, so the training sample
  # ends in 2017 and the forecast for 2018 is a one-step ahead forecast
  expect_length(external, 2)
  expect_equal(stats::tsp(external[[1]][["data"]][["train"]][["y"]])[2], 2017)
  expect_equal(as.numeric(external[[1]][["posterior"]][["forecast"]]), c(1, 2))
  expect_equal(stats::tsp(external[[2]][["data"]][["train"]][["y"]])[2], 2018)
  expect_equal(as.numeric(external[[2]][["posterior"]][["forecast"]]), c(3, NA))
})

test_that("the functions of the estimation workflow leave external forecasts unchanged", {
  external <- create_external_forecast(ext_forecasts(value = 1), ext_reference(),
                                       n_ahead = 2, data_lag = 1)

  expect_identical(add_priors(external, coef = list(v_i = 1)), external)
  expect_identical(add_initial_values(external), external)
  expect_identical(add_posterior_coefficients(external), external)
  expect_identical(add_posterior_loglik(external), external)
  expect_identical(add_forecast_input(external, n_ahead = 8), external)
  expect_identical(add_posterior_forecasts(external), external)
  expect_identical(thin(external, thin = 2), external)
})

test_that("selection_criteria only provides out-of-sample statistics", {
  full <- stats::window(diff(log(bvartools::e1)) * 100, end = c(1978, 4))
  external <- create_external_forecast(ext_forecasts(value = 0), ext_reference(),
                                       n_ahead = 2, data_lag = 1)
  external <- add_forecast_errors(external, test_sample = full)

  criteria <- selection_criteria(external)

  expect_s3_class(criteria, "selcrit")
  expect_named(criteria, c("model", "FE", "AFE", "RSFE"))
  expect_null(criteria[["LL"]])
  # A forecast of zero has an absolute forecast error equal to the realisation
  expect_true(all(criteria[["AFE"]][, "mean"] > 0))
})

test_that("external forecasts can be compared with models in one list", {
  full <- stats::window(diff(log(bvartools::e1)) * 100, end = c(1978, 4))
  external <- create_external_forecast(ext_forecasts(value = 0), ext_reference(),
                                       n_ahead = 2, data_lag = 1)

  models <- combine_models(add_posterior_loglik(ext_reference()), external)
  models <- add_forecast_errors(models, test_sample = full)
  criteria <- selection_criteria(models)

  expect_length(criteria, 2)
  # Only the model has in-sample criteria, which the printed comparison and the
  # choice of the best model have to tolerate
  expect_null(criteria[[2]][["LL"]])
  expect_output(print(criteria, relative = 1), "Model 2")
  expect_identical(choose_best_model(criteria, criterion = "BIC"), 1L)
})

test_that("a forecaster without a matching training sample is reported and dropped", {
  usable <- ext_forecasts(value = 1)
  usable[["forecaster"]] <- "A"
  # Publications long before the first training sample
  other <- usable
  other[["forecaster"]] <- "B"
  other[, "origin"] <- other[, "origin"] - 100

  expect_warning(external <- create_external_forecast(rbind(usable, other),
                                                      ext_reference(), n_ahead = 2,
                                                      by = "forecaster", data_lag = 1),
                 "'B'")
  expect_named(external, "A")
})

test_that("forecasts that cannot be matched to a training sample are reported", {
  forecasts <- ext_forecasts(value = 1)
  # Publications long before the first training sample
  forecasts[, "origin"] <- forecasts[, "origin"] - 100

  expect_error(create_external_forecast(forecasts, ext_reference(), n_ahead = 2),
               "could be matched")
})

test_that("input is validated", {
  expect_error(create_external_forecast(list(), ext_reference(), n_ahead = 2),
               "data frame")
  expect_error(create_external_forecast(ext_forecasts(), ext_reference(),
                                        n_ahead = 2, period = "date"),
               "'date'")
  expect_error(create_external_forecast(ext_forecasts(), ext_reference(),
                                        n_ahead = 2, select = "middle"),
               "'last' or 'first'")
})
