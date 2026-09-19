# Adding to a posterior returns a new object and leaves its argument alone.
#
# The C++ functions behind add_posterior_forecasts() and add_posterior_loglik()
# take the model as an Rcpp::List, which wraps the caller's R object rather than
# a copy of it. They used to write their result into that list, so the caller's
# object changed without an assignment. They also appended rather than set the
# element by name, so a second forecast left two `forecast` elements behind and
# `$forecast` kept returning the first, stale one.
#
# The object is compared through serialize() because a plain copy made in R
# would share its memory with the argument and change along with it.

test_that("add_posterior_forecasts leaves the object it is given unchanged", {
  object <- add_forecast_input(fx_var_fitted(), n_ahead = 2)
  before <- serialize(object, NULL)

  set.seed(1)
  result <- add_posterior_forecasts(object)

  expect_identical(serialize(object, NULL), before)
  expect_named(object[["posterior"]], c("a", "u_sigma_inv", "loglik"))
  expect_named(result[["posterior"]], c("a", "u_sigma_inv", "loglik", "forecast"))
})

test_that("add_posterior_loglik leaves the object it is given unchanged", {
  object <- fx_var_fitted()
  object[["posterior"]][["loglik"]] <- NULL
  before <- serialize(object, NULL)

  result <- add_posterior_loglik(object)

  expect_identical(serialize(object, NULL), before)
  expect_named(object[["posterior"]], c("a", "u_sigma_inv"))
  expect_named(result[["posterior"]], c("a", "u_sigma_inv", "loglik"))
})

test_that("add_posterior_loglik leaves an error correction model unchanged", {
  object <- fx_vec_fitted()
  object[["posterior"]][["loglik"]] <- NULL
  before <- serialize(object, NULL)

  result <- add_posterior_loglik(object)

  expect_identical(serialize(object, NULL), before)
  expect_false("loglik" %in% names(object[["posterior"]]))
  expect_true("loglik" %in% names(result[["posterior"]]))
})

test_that("forecasting an object that already carries a forecast replaces it", {
  object <- add_forecast_input(fx_var_fitted(), n_ahead = 2)

  set.seed(1)
  first <- add_posterior_forecasts(object)
  set.seed(2)
  second <- add_posterior_forecasts(first)
  set.seed(2)
  fresh <- add_posterior_forecasts(object)

  expect_identical(sum(names(second[["posterior"]]) == "forecast"), 1L)
  # The draws follow the seed: the second forecast is not the first one
  # returned again, and it is what the same seed produces from scratch.
  expect_false(identical(unclass(first[["posterior"]][["forecast"]][["forecasts"]]),
                         unclass(second[["posterior"]][["forecast"]][["forecasts"]])))
  expect_identical(second[["posterior"]][["forecast"]][["forecasts"]],
                   fresh[["posterior"]][["forecast"]][["forecasts"]])
})

test_that("adding the log likelihood again replaces the stored one", {
  object <- fx_var_fitted()

  result <- add_posterior_loglik(object)

  expect_identical(sum(names(result[["posterior"]]) == "loglik"), 1L)
  expect_equal(result[["posterior"]][["loglik"]], object[["posterior"]][["loglik"]])
})
