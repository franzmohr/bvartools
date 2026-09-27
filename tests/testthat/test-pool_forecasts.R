# Tests of pool_forecasts(): the equal-weight pool of the forecasts of several
# models. Checks that the pool mixes its members with equal weights, keeps the
# variables, horizons and forecast origins its members share, pools the log
# predictive densities as the density of the mixture, refuses what cannot be
# pooled, and passes through the workflow of a list of models. How forecasts
# are scored is left to test-selection_criteria.R and test-forecasts.R.

# A member whose forecasts are the constant 'value' in every draw, built from
# the expanding window fixture, so that the pooled values are known exactly.
pool_constant_member <- function(value, n_draws = NULL) {
  member <- fx_expanding_forecast()
  for (i in seq_along(member)) {
    x <- as.matrix(member[[i]][["posterior"]][["forecast"]][["forecasts"]])
    if (!is.null(n_draws)) {
      x <- x[seq_len(n_draws), , drop = FALSE]
    }
    x[] <- value
    member[[i]][["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(x)
    member[[i]][["posterior"]][["forecast"]][["errors"]] <- NULL
  }
  member
}

test_that("a pool of expanding windows is one pooled forecast per window", {
  set.seed(1)
  member <- fx_expanding_forecast()
  pool <- pool_forecasts(a = member, b = member)

  expect_s3_class(pool, "forecastpool")
  expect_s3_class(pool, "expandingwindow")
  expect_length(pool, length(member))
  expect_s3_class(pool[[1]], "poolwindow")
  expect_identical(pool[[1]][["model"]][["members"]], c("a", "b"))
  expect_identical(pool[[1]][["model"]][["endogen"]], member[[1]][["model"]][["endogen"]])

  # Forecasts made at the end of the same estimation sample are pooled
  ends <- function(x) vapply(x, function(w) stats::tsp(w[["data"]][["train"]][["y"]])[2], numeric(1))
  expect_equal(unname(ends(pool)), unname(ends(member)))
})

test_that("every member has the same weight in the pool", {
  set.seed(2)
  pool <- pool_forecasts(zero = pool_constant_member(0), ten = pool_constant_member(10))
  draws <- as.matrix(pool[[1]][["posterior"]][["forecast"]][["forecasts"]])

  # The mixture of a point mass at 0 and one at 10 with equal weights: half of
  # the draws are 0, the other half 10, and the mean is 5
  expect_equal(mean(draws == 0), 0.5)
  expect_equal(mean(draws), 5)
})

test_that("the weights stay equal when the members have different numbers of draws", {
  set.seed(3)
  pool <- pool_forecasts(pool_constant_member(0), pool_constant_member(10, n_draws = 4))
  draws <- as.matrix(pool[[1]][["posterior"]][["forecast"]][["forecasts"]])

  # Four draws from each, as many as the member with fewer draws has, not all
  # draws of the larger member, which would give it the larger weight
  expect_identical(nrow(draws), 8L)
  expect_identical(pool[[1]][["model"]][["draws"]], 4L)
  expect_equal(mean(draws), 5)
})

test_that("the pool keeps the variables and horizons all members share", {
  set.seed(4)
  member <- fx_expanding_forecast()
  smaller <- member
  keep <- c("Dp", "y")
  for (i in seq_along(smaller)) {
    endogen <- smaller[[i]][["model"]][["endogen"]]
    k <- length(endogen)
    x <- as.matrix(smaller[[i]][["posterior"]][["forecast"]][["forecasts"]])
    # One horizon only, and the variables in a different order
    cols <- match(keep, endogen)
    smaller[[i]][["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(x[, cols, drop = FALSE])
    smaller[[i]][["model"]][["endogen"]] <- keep
    smaller[[i]][["model"]][["k"]] <- length(keep)
  }
  pool <- pool_forecasts(member, smaller)

  expect_setequal(pool[[1]][["model"]][["endogen"]], keep)
  expect_identical(pool[[1]][["model"]][["h"]], 1)
  expect_identical(ncol(pool[[1]][["posterior"]][["forecast"]][["forecasts"]]), 2L)

  # Each column holds the draws of the same variable from both members
  big <- as.matrix(member[[1]][["posterior"]][["forecast"]][["forecasts"]])
  pooled <- as.matrix(pool[[1]][["posterior"]][["forecast"]][["forecasts"]])
  first <- pool[[1]][["model"]][["endogen"]][1]
  expect_true(all(pooled[, 1] %in% big[, match(first, member[[1]][["model"]][["endogen"]])]))
})

test_that("forecasts from estimation samples not every member has are left out", {
  set.seed(5)
  member <- fx_expanding_forecast()
  shorter <- member[-1]
  class(shorter) <- class(member)

  expect_message(pool <- pool_forecasts(member, shorter), "left out")
  expect_length(pool, length(member) - 1)
})

test_that("the log predictive densities are pooled as the density of the mixture", {
  set.seed(6)
  member <- fx_expanding_forecast()
  a <- b <- member
  for (i in seq_along(member)) {
    a[[i]][["predictive"]] <- list(loglik = log(rep(c(0.1, 0.3), 10)), period = 2000 + i)
    b[[i]][["predictive"]] <- list(loglik = log(rep(0.6, 20)), period = 2000 + i)
  }
  pool <- pool_forecasts(a, b)
  pooled <- pool[[1]][["predictive"]][["loglik"]]

  # The density of the equal-weight mixture is the average of the members'
  # densities, (0.2 + 0.6) / 2, not a density of the pooled draws
  expect_equal(log(mean(exp(pooled))), log(0.4))
  expect_equal(pool[[1]][["predictive"]][["period"]], 2001)

  # A member without predictive densities leaves the pool without them
  pool <- pool_forecasts(a, member)
  expect_null(pool[[1]][["predictive"]])
})

test_that("a pool is scored and compared like its members", {
  set.seed(7)
  full <- var_data()
  pool <- pool_forecasts(zero = pool_constant_member(0), ten = pool_constant_member(10))
  pool <- add_forecast_errors(pool, test_sample = full)

  # Its errors are the realised values minus the pooled draws
  realised <- as.matrix(stats::window(full, start = stats::tsp(pool[[1]][["data"]][["train"]][["y"]])[2] + 0.25))[1, ]
  errors <- as.matrix(pool[[1]][["posterior"]][["forecast"]][["errors"]])
  expect_equal(colMeans(errors)[1:3], unname(realised[pool[[1]][["model"]][["endogen"]]]) - 5,
               ignore_attr = TRUE)

  all_models <- combine_models(fx_expanding_forecast(), pool)
  expect_s3_class(all_models, "modellist")
  sc <- selection_criteria(all_models)
  expect_length(sc, 2)
  expect_false(is.null(sc[[2]][["FE"]]))
})

test_that("the steps that estimate a model leave a pool unchanged", {
  set.seed(8)
  pool <- pool_forecasts(fx_expanding_forecast(), fx_expanding_forecast())
  models <- combine_models(fx_expanding_forecast(), pool)

  expect_identical(add_posterior_forecasts(models)[[2]], pool)
  expect_identical(add_posterior_coefficients(models)[[2]], pool)
  expect_identical(add_predictive_loglik(pool), pool)
  expect_identical(thin(pool), pool)
  expect_error(write_to_hdf5(pool, folder = tempdir()), "not written to a file")
})

test_that("a pool is reproducible with set.seed()", {
  member <- fx_expanding_forecast()
  set.seed(9)
  first <- pool_forecasts(member, member)
  set.seed(9)
  second <- pool_forecasts(member, member)
  expect_identical(first, second)
})

test_that("what cannot be pooled is refused with a reason", {
  member <- fx_expanding_forecast()

  expect_error(pool_forecasts(member), "at least two members")

  # A single model and an expanding window
  expect_error(pool_forecasts(member, member[[1]]), "either expanding window exercises or single models")

  # A pool as a member of a pool
  pool <- pool_forecasts(member, member)
  expect_error(pool_forecasts(pool, member), "cannot be a member of another pool")

  # Members without a shared variable
  other <- member
  for (i in seq_along(other)) {
    other[[i]][["model"]][["endogen"]] <- c("a", "b", "c")
  }
  expect_error(pool_forecasts(member, other), "share no endogenous variable")

  # A member without forecasts
  unforecast <- fx_expanding_window()
  expect_error(pool_forecasts(unforecast, unforecast), "no forecasts|same estimation sample")
})

test_that("single models are pooled into one pooled forecast", {
  set.seed(10)
  member <- fx_expanding_forecast()
  pool <- pool_forecasts(member[[1]], member[[1]])

  expect_s3_class(pool, "poolwindow")
  expect_identical(nrow(as.matrix(pool[["posterior"]][["forecast"]][["forecasts"]])),
                   2L * nrow(as.matrix(member[[1]][["posterior"]][["forecast"]][["forecasts"]])))
})
