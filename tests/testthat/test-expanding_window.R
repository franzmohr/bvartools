test_that("use_expanding_window builds one model per forecast origin", {
  windows <- use_expanding_window(fx_var_model(), start = c(1976, 1))

  expect_s3_class(windows, "expandingwindow")
  expect_true(all(vapply(windows, inherits, logical(1), "bvarmodel")))

  # The first window ends just before `start`, and one window is added per
  # period up to the end of the sample.
  n_after_start <- nrow(stats::window(fx_var_model()[["data"]][["train"]][["y"]],
                                      start = c(1976, 1)))
  expect_length(windows, n_after_start + 1)
})

test_that("the estimation samples grow by one period at a time", {
  windows <- use_expanding_window(fx_var_model(), start = c(1976, 1))
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
  windows <- use_expanding_window(model, start = c(1978, 1))
  windows <- add_initial_values(windows)
  set.seed(19)
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
