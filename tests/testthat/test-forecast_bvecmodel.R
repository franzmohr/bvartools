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
  expect_s3_class(from_vec[["posterior"]][["forecast"]][["forecasts"]], "mcmc")
  expect_equal(unclass(from_vec[["posterior"]][["forecast"]][["forecasts"]]),
               unclass(from_var[["posterior"]][["forecast"]][["forecasts"]]), ignore_attr = TRUE)

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

    forecast <- simulated[["posterior"]][["forecast"]][["forecasts"]]
    expect_identical(dim(forecast), c(as.integer(fx_iterations), as.integer(4 * k)))
    expect_true(all(is.finite(forecast)))
    expect_identical(held[["model"]][["forecast_states"]], "hold")

    # From the same seed, the drift of the states is the one thing separating
    # the two.
    expect_false(isTRUE(all.equal(unclass(forecast),
                                  unclass(held[["posterior"]][["forecast"]][["forecasts"]]),
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
  set.seed(23)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)
  windows <- add_posterior_loglik(windows)

  windows <- add_forecast_input(windows, n_ahead = 2)
  windows <- add_posterior_forecasts(windows)
  expect_true(all(vapply(windows, function(x) !is.null(x[["posterior"]][["forecast"]][["forecasts"]]),
                         logical(1))))

  windows <- add_forecast_errors(windows, test_sample = full)
  expect_false(is.null(get_forecast_errors(windows[[1]])))

  criteria <- selection_criteria(windows)
  expect_true(all(c("FE", "AFE", "RSFE") %in% names(criteria)))
})

# A time varying VEC model on log levels times 100, whose error correction term
# is hundreds of times its typical change per period away from zero.
fx_vec_tvp_levels <- function(scaled = FALSE) {
  levels <- stats::window(at_domestic()[, c("y", "lr")], end = c(2005, 4)) * 100
  model <- create_bvecmodel(levels, p = 2, r = 1, const = "unrestricted", tvp = TRUE,
                            iterations = fx_iterations, burnin = fx_burnin)
  if (scaled) {
    model <- scale_error_correction(model, centre = TRUE)
  }
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                      coint = list(rho = 0.999), sigma = list(df = "k", scale = 1))
  model <- add_posterior_coefficients(add_seed(add_initial_values(model), 16180))
  if (scaled) {
    model <- rescale_error_correction(model)
  }
  add_forecast_input(model, n_ahead = 2)
}

test_that("simulating the cointegration vectors of series far from zero draws a warning", {
  vec <- fx_vec_tvp_levels()

  expect_warning(simulated <- add_posterior_forecasts(vec), "error correction term")
  expect_true(all(is.finite(simulated[["posterior"]][["forecast"]][["forecasts"]])))
  expect_no_warning(add_posterior_forecasts(vec, forecast_states = "hold"))

  # Interest rates and inflation in percent sit near zero and do not draw it.
  expect_no_warning(add_posterior_forecasts(add_forecast_input(fx_vec_tvp_fitted("gamma"), n_ahead = 2)))
})

test_that("a time varying VEC model estimated on centred series forecasts only with held states", {
  vec <- fx_vec_tvp_levels(scaled = TRUE)
  expect_true(vec[["model"]][["ect_rescaled"]])

  expect_error(add_posterior_forecasts(vec), "scale_error_correction")
  held <- add_posterior_forecasts(vec, forecast_states = "hold")
  expect_true(all(is.finite(held[["posterior"]][["forecast"]][["forecasts"]])))

  realised <- stats::window(at_domestic()[, c("y", "lr")], start = c(2006, 1), end = c(2006, 2)) * 100
  expect_no_error(add_predictive_loglik(held, test_sample = realised))
  simulate <- held
  simulate[["model"]][["forecast_states"]] <- "simulate"
  expect_error(add_predictive_loglik(simulate, test_sample = realised), "scale_error_correction")

  # The flag travels with the model through a file.
  path <- temp_h5_file()
  write_to_hdf5(vec, path)
  back <- read_model_from_hdf5(path)
  expect_error(add_posterior_forecasts(back), "scale_error_correction")
})

test_that("rescaling a constant VEC model leaves its forecasts free", {
  levels <- stats::window(at_domestic()[, c("y", "lr")], end = c(2005, 4)) * 100
  model <- create_bvecmodel(levels, p = 2, r = 1, const = "unrestricted",
                            iterations = fx_iterations, burnin = fx_burnin)
  model <- scale_error_correction(model, centre = TRUE)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0), coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1))
  model <- rescale_error_correction(add_posterior_coefficients(add_seed(add_initial_values(model), 16180)))
  expect_null(model[["model"]][["ect_rescaled"]])
  expect_no_warning(add_posterior_forecasts(add_forecast_input(model, n_ahead = 2)))
})
