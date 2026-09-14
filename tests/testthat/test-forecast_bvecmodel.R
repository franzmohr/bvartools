test_that("a constant VEC model forecasts in levels as its VAR representation does", {
  vec <- add_forecast_input(fx_vec_fitted(), n_ahead = 3)
  var <- add_forecast_input(vec_to_var(fx_vec_fitted()), n_ahead = 3)

  # The regressors of the forecast periods are the level VAR's.
  expect_equal(vec[["model"]][["h"]], var[["model"]][["h"]])
  expect_equal(vec[["data"]][["forecast"]][["x"]], var[["data"]][["forecast"]][["x"]])

  # With constant coefficients the two routes are the same simulation from the
  # same draws, so the same seed gives the same forecasts.
  set.seed(7)
  from_vec <- add_posterior_forecasts(vec)
  set.seed(7)
  from_var <- add_posterior_forecasts(var)

  expect_s3_class(from_vec, "bvecmodel")
  expect_s3_class(from_vec[["posterior"]][["forecast"]], "mcmc")
  expect_equal(unclass(from_vec[["posterior"]][["forecast"]]),
               unclass(from_var[["posterior"]][["forecast"]]), ignore_attr = TRUE)

  pred_vec <- stats::predict(from_vec)
  pred_var <- stats::predict(from_var)
  expect_s3_class(pred_vec, "bvarprd")
  expect_equal(pred_vec[["fcst"]], pred_var[["fcst"]])
  expect_equal(pred_vec[["y"]], pred_var[["y"]])
})

test_that("a time varying VEC model simulates its states forward unless told to hold them", {
  for (error in c("gamma", "sv")) {
    vec <- add_forecast_input(fx_vec_tvp_fitted(error), n_ahead = 4)
    k <- vec[["model"]][["k"]]

    set.seed(11)
    simulated <- add_posterior_forecasts(vec)
    set.seed(11)
    held <- add_posterior_forecasts(vec, forecast_states = "hold")

    forecast <- simulated[["posterior"]][["forecast"]]
    expect_identical(dim(forecast), c(as.integer(fx_iterations), as.integer(4 * k)))
    expect_true(all(is.finite(forecast)))
    expect_identical(held[["model"]][["forecast_states"]], "hold")

    # From the same seed, the drift of the states is the one thing separating
    # the two.
    expect_false(isTRUE(all.equal(unclass(forecast),
                                  unclass(held[["posterior"]][["forecast"]]),
                                  check.attributes = FALSE)))

    expect_s3_class(stats::predict(simulated), "bvarprd")
  }

  expect_error(add_posterior_forecasts(add_forecast_input(fx_vec_tvp_fitted("gamma"), n_ahead = 2),
                                       forecast_states = "sideways"))
})

test_that("a VEC model needs forecast input before it is forecast", {
  expect_error(add_posterior_forecasts(fx_vec_fitted()), "forecast horizon 'h'")
  expect_error(stats::predict(fx_vec_fitted()), "Missing element object\\$posterior\\$forecast")
})

test_that("VEC models over an expanding window are forecast and evaluated directly", {
  full <- vec_data()
  train <- stats::window(full, end = c(1994, 4))
  model <- create_bvecmodel(train, p = 2, r = 1, const = "unrestricted",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1))
  windows <- use_expanding_window(model, start = c(1994, 2))
  windows <- add_initial_values(windows)
  set.seed(23)
  windows <- add_posterior_coefficients(windows)
  windows <- add_posterior_loglik(windows)

  windows <- add_forecast_input(windows, n_ahead = 2)
  windows <- add_posterior_forecasts(windows)
  expect_true(all(vapply(windows, function(x) !is.null(x[["posterior"]][["forecast"]]),
                         logical(1))))

  windows <- add_forecast_errors(windows, test_sample = full)
  expect_false(is.null(get_forecast_errors(windows[[1]])))

  criteria <- selection_criteria(windows)
  expect_true(all(c("FE", "AFE", "RSFE") %in% names(criteria)))
})
