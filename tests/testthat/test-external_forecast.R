# The reference model of the comparisons below. The training samples of its
# windows end in 1996Q1 to 1997Q1, the test data run to 1998Q1.
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
  full <- var_data()
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
  forecast <- external[[1]][["posterior"]][["forecast"]][["forecasts"]]

  # One row, because external forecasts are point forecasts, and one column per
  # endogenous variable within each of the two forecast horizons
  expect_identical(dim(forecast), c(1L, 6L))
  expect_true(all(forecast == 1))
})

test_that("perfect foresight leaves no forecast errors", {
  full <- var_data()
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

  expect_identical(ncol(external[[1]][["posterior"]][["forecast"]][["forecasts"]]), 6L)
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
  expect_true(all(last[[1]][["posterior"]][["forecast"]][["forecasts"]] == 2))
  expect_true(all(first[[1]][["posterior"]][["forecast"]][["forecasts"]] == 1))
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
  expect_equal(as.numeric(external[[1]][["posterior"]][["forecast"]][["forecasts"]]), c(1, 2))
  expect_equal(stats::tsp(external[[2]][["data"]][["train"]][["y"]])[2], 2018)
  expect_equal(as.numeric(external[[2]][["posterior"]][["forecast"]][["forecasts"]]), c(3, NA))
})

test_that("a publication on the first day of a quarter belongs to that quarter", {
  quarterly <- stats::ts(matrix(seq_len(40) / 10, 40, 1, dimnames = list(NULL, "gdp")),
                         start = c(2010, 1), frequency = 4)
  model <- create_bvarmodel(quarterly, p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- use_expanding_window(model, start = c(2018, 1))

  # 1 April is day 91 of 365. As a share of the year it fell short of 0.25 and
  # counted as the first quarter, so with one quarter of data lag it was matched
  # to the window ending in the fourth quarter of the year before, like 31 March.
  forecasts <- data.frame(origin = as.Date(c("2019-03-31", "2019-04-01")),
                          period = c(2019.25, 2019.5), variable = "gdp", value = c(1, 2))
  external <- create_external_forecast(forecasts, model, n_ahead = 2, data_lag = 1)

  expect_length(external, 2)
  expect_equal(stats::tsp(external[[1]][["data"]][["train"]][["y"]])[2], 2018.75)
  expect_equal(stats::tsp(external[[2]][["data"]][["train"]][["y"]])[2], 2019)

  # As character strings a date is read as the first of its month, which put
  # every publication in April, July or October in the quarter before.
  forecasts <- data.frame(origin = c("2019-07-20", "2019-10-20"),
                          period = c(2019.5, 2019.75), variable = "gdp", value = c(1, 2))
  external <- create_external_forecast(forecasts, model, n_ahead = 2, data_lag = 1)
  expect_equal(unname(sapply(external, function(x) stats::tsp(x[["data"]][["train"]][["y"]])[2])),
               c(2019.25, 2019.5))
})

test_that("forecasts at another frequency than the data are refused", {
  # The reference model is quarterly. Annual periods used to be rounded to the
  # first quarter of each year, and an annual growth rate scored as a quarterly
  # one.
  annual <- expand.grid(origin = c(1996.3, 1996.55, 1996.8, 1997.05, 1997.3),
                        period = 1997:1998, variable = dimnames(var_data())[[2]],
                        stringsAsFactors = FALSE)
  annual[["value"]] <- 1
  expect_error(create_external_forecast(annual, ext_reference(), n_ahead = 8),
               "appear to be annual forecasts")

  # As calendar dates, which are rounded rather than read as numbers
  annual_dates <- annual
  annual_dates[["period"]] <- paste0(annual_dates[["period"]], "-12-31")
  expect_error(create_external_forecast(annual_dates, ext_reference(), n_ahead = 8),
               "appear to be annual forecasts")

  # One period per publication gives no spacing to read, so the position within
  # the year does: every period is a first quarter.
  single <- annual[annual[["period"]] == floor(annual[["origin"]]) + 1, ]
  expect_error(create_external_forecast(single, ext_reference(), n_ahead = 8),
               "appear to be annual forecasts")

  # Monthly forecasts, three of which fall into the same quarter
  monthly <- expand.grid(origin = 1996.3, period = 1996.5 + (0:5) / 12,
                         variable = "y", stringsAsFactors = FALSE)
  monthly[["value"]] <- 1
  expect_error(create_external_forecast(monthly, ext_reference(), n_ahead = 4),
               "same period of the data, which is quarterly")

  # Quarterly forecasts of one period each, one publication per quarter, are
  # at the frequency of the data.
  quarterly <- ext_forecasts(value = 1)
  quarterly <- quarterly[abs(quarterly[["period"]] - quarterly[["origin"]] + 0.05) < 1e-8, ]
  expect_equal(nrow(quarterly), length(ext_ends()) * ncol(var_data()))
  expect_s3_class(create_external_forecast(quarterly, ext_reference(), n_ahead = 2),
                  "externalforecast")
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
  full <- var_data()
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
  full <- var_data()
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

# --- annual forecasts of a quarterly model ------------------------------------

# An expanding window of the quarterly reference, whose eight quarters cover the
# year of every forecast origin and the next. The first training sample ends in
# 1995Q1 and the last, like the data, in 1997Q1, so the realised annual figures
# run to 1996.
ann_windows <- function() {
  cached_fixture("annual_windows", {
    train <- stats::window(var_data(), end = c(1997, 1))
    model <- create_bvarmodel(train, p = 1, deterministic = "const",
                              iterations = 10, burnin = 5)
    model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                        sigma = list(df = 1, scale = 0.0001))
    windows <- use_expanding_window(model, start = c(1995, 2))
    windows <- add_initial_values(windows)
    windows <- add_posterior_coefficients(add_seed(windows, 20200))
    windows <- add_forecast_input(windows, n_ahead = 8)
    add_posterior_forecasts(windows)
  })
}

# Output growth in percent, the log change of a level times 100, and inflation
# as it is
ann_code <- c(y = 5, Dp = 1)
ann_kind <- c(y = "growth", Dp = "level")
ann_aggregate <- function(object, ...) {
  aggregate_forecasts(object, code = ann_code, scale = 100, ...)
}

# The annual figure of a year from quarterly values in the units of the model:
# the growth of the annual average of the levels, whose log levels in percent
# are the cumulated changes, or the average over the year.
ann_by_hand <- function(time, values, year, type) {
  year_of <- floor(time + 1e-8)
  if (type == "level") {
    return(mean(values[year_of == year]))
  }
  levels <- exp(cumsum(values) / 100)
  100 * (mean(levels[year_of == year]) / mean(levels[year_of == year - 1]) - 1)
}

test_that("each draw of a quarterly forecast becomes a draw of the annual figures", {
  windows <- ann_windows()
  annual <- ann_aggregate(windows)

  expect_s3_class(annual, "expandingwindow")
  expect_length(annual, length(windows))

  for (i in c(1, 3, length(windows))) {
    window_i <- windows[[i]]
    end <- stats::tsp(window_i[["data"]][["train"]][["y"]])[2]
    observed <- stats::window(window_i[["data"]][["original"]][["endogen"]], end = end)
    draws <- window_i[["posterior"]][["forecast"]][["forecasts"]]
    result <- annual[[i]][["posterior"]][["forecast"]][["forecasts"]]

    expect_s3_class(result, "mcmc")
    expect_identical(dimnames(result)[[2]], c("y_1", "Dp_1", "y_2", "Dp_2"))
    expect_identical(nrow(result), nrow(draws))

    # Horizon 1 is the year of the first quarter after the training sample,
    # whose quarters up to the end of the sample are observed
    first_year <- floor(end + 0.25 + 1e-8)
    time <- c(stats::time(observed), end + (1:8) / 4)
    for (s in c(1, nrow(draws))) {
      for (v in names(ann_code)) {
        values <- c(observed[, v], draws[s, (0:7) * 3 + match(v, dimnames(observed)[[2]])])
        for (h in 1:2) {
          expect_equal(unname(result[s, paste0(v, "_", h)]),
                       ann_by_hand(time, values, first_year + h - 1, ann_kind[[v]]))
        }
      }
    }
  }
})

# Seven untransformed quarterly series from 1990 to 2001, one per transformation
# code of FRED-QD: the positive levels of a growing series for the
# multiplicative codes and a wandering rate for the additive ones.
ann_levels <- function() {
  set.seed(1990)
  n <- 48
  growing <- function() 100 * exp(cumsum(0.005 + 0.01 * stats::rnorm(n)))
  wandering <- function() 5 + cumsum(0.2 * stats::rnorm(n))
  x <- cbind(c1 = wandering(), c2 = wandering(), c3 = wandering(), c4 = growing(),
             c5 = growing(), c6 = growing(), c7 = growing())
  stats::ts(x, start = 1990, frequency = 4)
}
ann_codes <- c(c1 = 1, c2 = 2, c3 = 3, c4 = 4, c5 = 5, c6 = 6, c7 = 7)

# An expanding window over the transformed series, whose training samples end in
# 1997Q4 to 1999Q4, and whose forecast draws are all replaced by the transformed
# values that followed, eight quarters of them. The data of the models end in
# 1999Q4, the levels in 2001Q4.
ann_foresight <- function(scale = 1) {
  data <- transform_variables(ann_levels(), ann_codes)
  data[, 4:7] <- data[, 4:7] * scale
  data <- stats::window(data, start = c(1990, 3))
  model <- create_bvarmodel(stats::window(data, end = c(1999, 4)), p = 1,
                            deterministic = "const", iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0.1, v_i_det = 0.01),
                      sigma = list(df = "k", scale = 1))
  windows <- use_expanding_window(model, start = c(1998, 1))
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(add_seed(windows, 1990))
  windows <- add_forecast_input(windows, n_ahead = 8)
  windows <- add_posterior_forecasts(windows)
  for (i in seq_along(windows)) {
    end <- stats::tsp(windows[[i]][["data"]][["train"]][["y"]])[2]
    future <- stats::window(data, start = end + 0.25, end = end + 2)
    draws <- windows[[i]][["posterior"]][["forecast"]][["forecasts"]]
    draws[] <- rep(as.numeric(t(future)), each = nrow(draws))
    windows[[i]][["posterior"]][["forecast"]][["forecasts"]] <- draws
  }
  windows
}

test_that("perfect foresight of the transformed series is perfect foresight of the annual figures", {
  for (scale in c(1, 100)) {
    windows <- ann_foresight(scale)
    for (target in c("average", "q4q4")) {
      # Reversing the transformation of every code and aggregating the result
      # recovers the annual figures of the untransformed series
      annual <- aggregate_forecasts(windows, code = ann_codes, target = target,
                                    levels = ann_levels(), scale = scale)
      annual <- add_forecast_errors(annual)
      errors <- unlist(lapply(annual, get_forecast_errors))
      # Every window forecasts two years, which the levels cover
      expect_true(all(sapply(annual, function(x) ncol(get_forecast_errors(x))) ==
                        2 * length(ann_codes)))
      expect_equal(errors, rep(0, length(errors)), ignore_attr = TRUE)

      # Codes 1, 4 and 5 are reversed from the data of the models alone
      from_data <- aggregate_forecasts(windows, code = ann_codes[c(1, 4, 5)],
                                       target = target, scale = scale)
      with_levels <- aggregate_forecasts(windows, code = ann_codes[c(1, 4, 5)],
                                         target = target, levels = ann_levels(),
                                         scale = scale)
      for (i in seq_along(windows)) {
        expect_equal(from_data[[i]][["posterior"]][["forecast"]][["forecasts"]],
                     with_levels[[i]][["posterior"]][["forecast"]][["forecasts"]])
      }
      # Their realised values end with the data of the models, which is where
      # the levels reach further
      expect_equal(stats::tsp(from_data[[1]][["data"]][["original"]][["endogen"]])[2], 1999)
      expect_equal(stats::tsp(with_levels[[1]][["data"]][["original"]][["endogen"]])[2], 2001)
    }
  }
})

test_that("the two targets are the annual average and the fourth quarter", {
  x <- ann_levels()
  in_year <- function(v, year) as.numeric(stats::window(x[, v], start = year, end = year + 0.75))
  average <- aggregate_forecasts(ann_foresight(), code = ann_codes, levels = x)
  q4q4 <- aggregate_forecasts(ann_foresight(), code = ann_codes, levels = x,
                              target = "q4q4")

  # The first window ends in 1997Q4, so horizon 2 is 1999
  average <- average[[1]][["posterior"]][["forecast"]][["forecasts"]][1, ]
  q4q4 <- q4q4[[1]][["posterior"]][["forecast"]][["forecasts"]][1, ]
  expect_equal(unname(average["c2_2"]), mean(in_year("c2", 1999)))
  expect_equal(unname(average["c6_2"]),
               100 * (mean(in_year("c6", 1999)) / mean(in_year("c6", 1998)) - 1))
  expect_equal(unname(q4q4["c2_2"]), in_year("c2", 1999)[4])
  expect_equal(unname(q4q4["c6_2"]), 100 * (in_year("c6", 1999)[4] / in_year("c6", 1998)[4] - 1))
})

test_that("levels must reproduce the data of the models", {
  windows <- ann_foresight(scale = 100)
  # The data are in percent, the default scale is 1
  expect_error(aggregate_forecasts(windows, code = ann_codes, levels = ann_levels()),
               "does not reproduce the data of variable 'c4' under code 4 and scale 1")
  expect_error(aggregate_forecasts(windows, code = ann_codes, scale = 100),
               "'c2', 'c3', 'c6', 'c7'")
  expect_error(aggregate_forecasts(windows, code = ann_codes, scale = 100,
                                   levels = ann_levels()[, 1:6]),
               "'c7'")
  expect_error(aggregate_forecasts(windows, code = ann_codes, scale = 100,
                                   levels = as.data.frame(ann_levels())),
               "class 'ts'")
})

test_that("aggregated models carry the realised annual figures they are scored against", {
  annual <- ann_aggregate(ann_windows())
  full <- var_data()

  # The first window ends in 1995Q1: horizon 1 is 1995, horizon 2 is 1996
  realised <- annual[[1]][["data"]][["test"]][["y"]]
  time <- as.numeric(stats::time(full))
  expect_equal(unname(realised[, "y"]),
               sapply(1995:1996, function(year) ann_by_hand(time, full[, "y"], year, "growth")))
  expect_equal(unname(realised[, "Dp"]),
               sapply(1995:1996, function(year) ann_by_hand(time, full[, "Dp"], year, "level")))

  # The data end in 1997Q1, so the windows of 1996Q4 and 1997Q1 forecast years
  # without a realised figure, and are left without forecast errors
  last <- annual[[length(annual)]]
  expect_null(last[["data"]][["test"]][["y"]])
  expect_identical(nrow(annual[[length(annual) - 2]][["data"]][["test"]][["y"]]), 1L)

  scored <- add_forecast_errors(annual)
  errors <- get_forecast_errors(scored[[1]])
  forecasts <- scored[[1]][["posterior"]][["forecast"]][["forecasts"]]
  expect_equal(unname(errors[1, ]), unname(as.numeric(t(realised)) - forecasts[1, ]))

  expect_null(get_forecast_errors(scored[[length(scored)]]))

  criteria <- selection_criteria(scored)
  expect_identical(criteria[["RSFE"]][["variable"]], c("y", "Dp", "y", "Dp"))
  expect_identical(criteria[["RSFE"]][["h"]], c(1L, 1L, 2L, 2L))

  # A quarterly test sample would be matched by time alone, 1996 being a year
  # and the first quarter of 1996 alike
  expect_error(add_forecast_errors(annual, test_sample = full), "quarterly data")
})

test_that("annual external forecasts are matched to the quarterly training samples", {
  annual <- ann_aggregate(ann_windows())

  # Publications in the second and the fourth quarter of 1995 for the current
  # and the next year, which see the data up to the quarter before
  forecasts <- expand.grid(origin = c(1995.3, 1995.8), year = 0:1,
                           variable = names(ann_code), stringsAsFactors = FALSE)
  forecasts[["period"]] <- 1995 + forecasts[["year"]]
  forecasts[["value"]] <- forecasts[["origin"]] + forecasts[["year"]]

  external <- create_external_forecast(forecasts, annual, data_lag = 1)
  expect_length(external, 2)
  expect_identical(external[[1]][["model"]][["h"]], 2L)
  expect_equal(sapply(external, function(x) x[["model"]][["aggregation"]][["end"]]),
               c(1995, 1995.5))

  # Whichever quarter a forecast was published in, the current year is horizon 1
  for (i in 1:2) {
    values <- external[[i]][["posterior"]][["forecast"]][["forecasts"]]
    expect_equal(unname(values[1, ]), rep(c(1995.3, 1995.8)[i] + 0:1, each = 2))
    expect_equal(stats::tsp(external[[i]][["data"]][["train"]][["y"]])[2:3], c(1994, 1))
  }

  expect_output(print(external), "Annual forecasts of models estimated on quarterly data")
  expect_output(print(external), "1995.5")
})

test_that("annual external forecasts are scored against the annual figures of the models", {
  annual <- ann_aggregate(ann_windows())
  realised <- annual[[1]][["data"]][["original"]][["endogen"]]

  forecasts <- expand.grid(origin = 1995.3 + (0:8) / 4, year = 0:1,
                           variable = names(ann_code), stringsAsFactors = FALSE)
  forecasts[["period"]] <- floor(forecasts[["origin"]]) + forecasts[["year"]]
  forecasts[["value"]] <- mapply(function(period, variable) {
    value <- realised[round(stats::time(realised)) == period, variable]
    if (length(value) == 0) NA_real_ else value
  }, forecasts[["period"]], forecasts[["variable"]])

  # Perfect foresight of the annual figures leaves no forecast errors, although
  # the external forecasts come without a test sample of their own
  external <- create_external_forecast(forecasts, annual, data_lag = 1)
  external <- add_forecast_errors(external)
  errors <- unlist(lapply(external, get_forecast_errors))
  expect_true(length(errors) > 0)
  expect_equal(errors[!is.na(errors)], rep(0, sum(!is.na(errors))), ignore_attr = TRUE)

  # Models and forecasters are compared on the same annual basis
  models <- add_forecast_errors(combine_models(annual, external))
  criteria <- selection_criteria(models)
  expect_length(criteria, 2)
  expect_identical(criteria[[1]][["RSFE"]][, c("variable", "h")],
                   criteria[[2]][["RSFE"]][, c("variable", "h")])
  expect_equal(criteria[[2]][["RSFE"]][["mean"]], rep(0, 4))
})

test_that("the forecasts of a VEC model are aggregated from its levels", {
  vec <- add_posterior_forecasts(add_forecast_input(fx_vec_fitted(), n_ahead = 4))
  annual <- aggregate_forecasts(vec, code = c(lr = 1))

  expect_s3_class(annual, "bvarmodel")
  expect_false(inherits(annual, "bvecmodel"))

  levels <- vec_to_var(fx_vec_fitted())[["data"]][["train"]][["y"]]
  end <- stats::tsp(levels)[2]
  draws <- vec[["posterior"]][["forecast"]][["forecasts"]]
  year <- floor(end + 0.25 + 1e-8)
  observed <- levels[floor(stats::time(levels) + 1e-8) == year, "lr"]
  n_future <- 4 - length(observed)
  expect_equal(unname(annual[["posterior"]][["forecast"]][["forecasts"]][1, 1]),
               mean(c(observed, draws[1, (seq_len(n_future) - 1) * 2 + 1])))
})

test_that("annual aggregation refuses what it cannot aggregate", {
  windows <- ann_windows()

  expect_error(aggregate_forecasts(windows, code = c(5, 1)), "named")
  expect_error(aggregate_forecasts(windows, code = c(y = "growth")), "named numeric")
  expect_error(aggregate_forecasts(windows, code = c(y = 8)), "1 to 7")
  expect_error(aggregate_forecasts(windows, code = c(gdp = 5)), "'gdp'")
  expect_error(ann_aggregate(windows, target = "q1q1"), "'average' or 'q4q4'")
  expect_error(aggregate_forecasts(windows, code = ann_code, scale = 0), "positive")
  expect_error(ann_aggregate(fx_expanding_forecast()), "does not cover a whole year")
  expect_error(ann_aggregate(fx_var_fitted()), "add_posterior_forecasts")

  annual <- ann_aggregate(windows)
  expect_error(ann_aggregate(annual), "already aggregated")

  external <- create_external_forecast(ext_forecasts(value = 1), ext_reference(),
                                       n_ahead = 2, data_lag = 1)
  expect_error(ann_aggregate(external), "external forecasts")

  forecasts <- data.frame(origin = 1995.3, period = 1995, variable = "y", value = 1)
  expect_error(create_external_forecast(forecasts, combine_models(annual, windows)),
               "aggregate_forecasts")

  # Annual forecasts given to the quarterly models name the remedy
  annual_forecasts <- expand.grid(origin = 1995.3, period = 1995:1996,
                                  variable = "y", value = 1)
  expect_error(create_external_forecast(annual_forecasts, windows), "aggregate_forecasts")
})
