test_that("use_expanding_window builds one model per forecast origin", {
  windows <- use_expanding_window(fx_var_model(), start = c(1995, 2))

  expect_s3_class(windows, "expandingwindow")
  expect_true(all(vapply(windows, inherits, logical(1), "bvarmodel")))

  # The first window ends just before `start`, and one window is added per
  # period up to the end of the sample.
  n_after_start <- nrow(stats::window(fx_var_model()[["data"]][["train"]][["y"]],
                                      start = c(1995, 2)))
  expect_length(windows, n_after_start + 1)
})

test_that("the estimation samples grow by one period at a time", {
  windows <- use_expanding_window(fx_var_model(), start = c(1995, 2))
  nobs <- vapply(windows,
                 function(x) nrow(x[["data"]][["train"]][["y"]]), integer(1))
  starts <- vapply(windows,
                   function(x) stats::tsp(x[["data"]][["train"]][["y"]])[1],
                   numeric(1))

  expect_equal(diff(nobs), rep(1L, length(nobs) - 1))
  # Expanding, not rolling: every window starts at the same period.
  expect_length(unique(starts), 1)
  # The last window is the full sample.
  expect_identical(nobs[length(nobs)],
                   nrow(fx_var_model()[["data"]][["train"]][["y"]]))
})

test_that("the whole workflow runs over an expanding window", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  windows <- use_expanding_window(model, start = c(1997, 2))
  set.seed(19)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)
  windows <- add_posterior_loglik(windows)

  expect_s3_class(windows, "expandingwindow")
  expect_true(all(vapply(
    windows, function(x) !is.null(x[["posterior"]][["a"]]), logical(1)
  )))
  expect_true(all(vapply(
    windows, function(x) nrow(x[["posterior"]][["a"]][["coeffs"]]), integer(1)
  ) == 10L))
})

test_that("expanding window results can be summarised and thinned", {
  windows <- fx_expanding_window()

  expect_s3_class(summary(windows), "list")
  thinned <- thin(windows, thin = 2)
  expect_s3_class(thinned, "expandingwindow")
  expect_true(all(vapply(
    thinned, function(x) nrow(x[["posterior"]][["a"]][["coeffs"]]), integer(1)
  ) == 5L))
})

test_that("forecasts can be produced for every window", {
  windows <- fx_expanding_forecast()
  forecasts <- stats::predict(windows, n_ahead = 2)

  expect_length(forecasts, length(windows))
  expect_true(all(vapply(forecasts, inherits, logical(1), "bvarprd")))
})

test_that("the out-of-sample criteria summarise the forecast errors", {
  criteria <- selection_criteria(fx_expanding_forecast())
  spec <- fx_expanding_forecast()[[1]][["model"]]

  expect_s3_class(criteria, "selcrit")
  expect_true(all(c("FE", "AFE", "RSFE") %in% names(criteria)))

  # One row per variable and forecast horizon.
  expect_identical(nrow(criteria[["FE"]]),
                   as.integer(spec[["k"]] * spec[["h"]]))
  expect_identical(criteria[["FE"]][["variable"]],
                   rep(spec[["endogen"]], spec[["h"]]))
  expect_identical(criteria[["AFE"]][["h"]],
                   rep(seq_len(spec[["h"]]), each = spec[["k"]]))

  # Absolute errors cannot be negative and cannot be smaller than the raw ones.
  expect_true(all(criteria[["AFE"]][["mean"]] >= 0))
  expect_true(all(criteria[["AFE"]][["mean"]] >= abs(criteria[["FE"]][["mean"]])))
})

test_that("the in-sample criteria of an expanding window need no forecasts", {
  windows <- fx_expanding_window()
  criteria <- selection_criteria(windows)

  expect_s3_class(criteria, "selcrit")
  expect_false(any(c("FE", "AFE", "RSFE") %in% names(criteria)))
  # They are the criteria of the last window, estimated on the most data.
  expect_equal(criteria[["LL"]], selection_criteria(windows[[length(windows)]])[["LL"]])
  expect_output(print(criteria), "In-sample")
})

test_that("an expanding window with neither log-likelihood nor forecast errors is refused", {
  windows <- fx_expanding_window()
  for (i in seq_along(windows)) {
    windows[[i]][["posterior"]][["loglik"]] <- NULL
  }
  expect_error(selection_criteria(windows), "log-likelihood")
})

test_that("the forecasts of every window start from the last period of that window", {
  # A structural VARX of the first differences of at_macrodata, a window() of
  # it, and a VEC model in levels through its VAR representation. The lags of
  # the first forecast period used to come from the end of the whole series,
  # whatever the estimation sample of the model.
  domestic <- diff(at_data()) * 100
  foreign <- diff(bvartools::at_macrodata[["foreign"]][, c("y.s", "poil")]) * 100
  svarx <- create_bvarmodel(stats::window(domestic, end = c(2019, 4)), p = 2,
                            exogen = stats::window(foreign, end = c(2019, 4)), s = 1,
                            deterministic = "const", structural = TRUE, error = "gamma",
                            iterations = 10, burnin = 5)

  windows <- add_forecast_input(use_expanding_window(svarx, start = c(2019, 1)),
                                n_ahead = 2, exogen = foreign)
  expect_length(windows, 5)
  for (w in windows) {
    y <- as.matrix(w[["data"]][["train"]][["y"]])
    n <- nrow(y)
    x <- w[["data"]][["forecast"]][["x"]]
    expect_equal(as.numeric(x[1, 1:6]), as.numeric(c(y[n, ], y[n - 1, ])))
    expect_equal(as.numeric(x[2, 4:6]), as.numeric(y[n, ]))
    ends <- stats::tsp(w[["data"]][["train"]][["y"]])[2]
    next_period <- which(abs(stats::time(foreign) - (ends + 0.25)) < 1e-6)
    expect_equal(as.numeric(x[1, 7:8]), as.numeric(foreign[next_period, ]))
  }

  cut <- add_forecast_input(window(svarx, end = c(2015, 4)), n_ahead = 1, exogen = foreign)
  y <- as.matrix(cut[["data"]][["train"]][["y"]])
  expect_equal(as.numeric(cut[["data"]][["forecast"]][["x"]][1, 1:3]), as.numeric(y[nrow(y), ]))

  vec <- create_bvecmodel(stats::window(at_data(), end = c(2019, 4)) * 100, p = 2, r = 1,
                          const = "unrestricted", iterations = 10, burnin = 5)
  vec <- add_priors(vec, coef = list(v_i = 0, v_i_det = 0), coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1e-4))
  set.seed(1)
  vec_windows <- add_initial_values(use_expanding_window(vec, start = c(2019, 3)))
  vec_windows <- add_forecast_input(vec_to_var(add_posterior_coefficients(vec_windows)), n_ahead = 1)
  for (w in vec_windows) {
    y <- as.matrix(w[["data"]][["train"]][["y"]])
    expect_equal(as.numeric(w[["data"]][["forecast"]][["x"]][1, 1:3]), as.numeric(y[nrow(y), ]))
  }
})
