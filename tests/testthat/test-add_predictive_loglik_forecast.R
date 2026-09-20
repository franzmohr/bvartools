# Scoring a forecast against what its horizon realised.
#
# Each column of posterior$forecast$loglik is the density of that period's
# realised observation given the ones before it, so the columns sum to the log
# predictive likelihood of the whole stretch rather than being h marginals. What
# the tests here pin is the shape it comes back in, that the realised history
# really reaches the later columns, that a VEC model is scored in levels, and
# that the numbers are the normal density they are claimed to be.

# The series the fixtures are built on, past the period they stop at, so that
# there is something for the forecast to be scored against.
var_full_data <- function() {
  levels <- at_domestic()
  stats::ts.intersect(y = diff(levels[, "y"]), Dp = levels[, "Dp"],
                      r = diff(levels[, "r"])) * 100
}

vec_full_data <- function() {
  at_domestic()[, c("lr", "Dp")] * 100
}

# A fitted fixture with a forecast of `n_ahead` periods.
with_forecast <- function(object, n_ahead = 4) {
  add_posterior_forecasts(add_forecast_input(object, n_ahead = n_ahead))
}


test_that("a forecast is scored against the periods its horizon realised", {
  object <- with_forecast(fx_var_fitted())
  scored <- add_predictive_loglik(object, test_sample = var_full_data())
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_s3_class(loglik, "mcmc")
  # One row per draw and one column per scored period, as every other block of
  # the posterior is laid out.
  expect_identical(dim(loglik), c(fx_iterations, 4L))
  expect_true(all(is.finite(loglik)))
  expect_identical(coda::mcpar(loglik),
                   coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]))

  # The paths are still there beside it: the two are members of one group.
  expect_true(all(c("forecasts", "loglik") %in%
                    names(scored[["posterior"]][["forecast"]])))
})


test_that("the score is the normal density of the realised values", {
  object <- with_forecast(fx_var_fitted())
  scored <- add_predictive_loglik(object, test_sample = var_full_data())

  k <- scored[["model"]][["k"]]
  p <- scored[["model"]][["p"]]
  y <- as.matrix(scored[["data"]][["test"]][["y"]])
  h <- nrow(y)

  # The regressors of the scored periods, with the lags taken from what was
  # realised rather than from what the forecast simulated. That substitution is
  # the whole content of the claim being checked.
  x <- as.matrix(scored[["data"]][["forecast"]][["x"]])[1:h, , drop = FALSE]
  for (i in 2:h) {
    for (j in 1:min(p, i - 1)) {
      x[i, ((j - 1) * k + 1):(j * k)] <- y[i - j, ]
    }
  }

  a <- as.matrix(scored[["posterior"]][["a"]][["coeffs"]])
  s <- as.matrix(scored[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  expected <- matrix(NA_real_, nrow(a), h)
  for (draw in seq_len(nrow(a))) {
    coeffs <- matrix(a[draw, ], k)
    precision <- matrix(s[draw, ], k, k)
    for (i in seq_len(h)) {
      u <- y[i, ] - coeffs %*% x[i, ]
      expected[draw, i] <- -k * log(2 * pi) / 2 +
        sum(log(diag(chol(precision)))) -
        as.numeric(t(u) %*% precision %*% u) / 2
    }
  }

  expect_equal(unclass(scored[["posterior"]][["forecast"]][["loglik"]]),
               expected, ignore_attr = TRUE)
})


test_that("a scored model carries what it was scored against", {
  object <- with_forecast(fx_var_fitted())
  scored <- add_predictive_loglik(object, test_sample = var_full_data())
  realised <- scored[["data"]][["test"]][["y"]]

  expect_identical(dim(realised), c(4L, as.integer(object[["model"]][["k"]])))

  # And is scored again from them, to the same numbers, without the sample.
  again <- add_predictive_loglik(scored)
  expect_equal(unclass(again[["posterior"]][["forecast"]][["loglik"]]),
               unclass(scored[["posterior"]][["forecast"]][["loglik"]]))
})


test_that("the realised history reaches the columns its lags belong to", {
  object <- with_forecast(fx_var_fitted())
  test_sample <- var_full_data()
  first <- add_predictive_loglik(object, test_sample = test_sample)

  moved <- test_sample
  starts_at <- stats::tsp(object[["data"]][["train"]][["y"]])[2] +
    1 / stats::tsp(object[["data"]][["train"]][["y"]])[3]
  moved[which(abs(stats::time(moved) - starts_at) < 1e-8), ] <-
    moved[which(abs(stats::time(moved) - starts_at) < 1e-8), ] + 5
  second <- add_predictive_loglik(object, test_sample = moved)

  a <- unclass(first[["posterior"]][["forecast"]][["loglik"]])
  b <- unclass(second[["posterior"]][["forecast"]][["loglik"]])

  # The period that moved, and the one after it, which conditions on it as a
  # lag. The fixture has one lag, so the third is untouched -- which is what
  # says the history enters as a lag rather than the whole column being
  # recomputed from a different object.
  expect_false(isTRUE(all.equal(a[, 1], b[, 1])))
  expect_false(isTRUE(all.equal(a[, 2], b[, 2])))
  expect_equal(a[, 3], b[, 3])
})


test_that("fewer realised periods than the horizon are scored on their own", {
  object <- with_forecast(fx_var_fitted(), n_ahead = 4)
  starts_at <- stats::tsp(object[["data"]][["train"]][["y"]])[2]
  short <- stats::window(var_full_data(), end = starts_at + 2 / 4)

  scored <- add_predictive_loglik(object, test_sample = short)
  expect_identical(ncol(scored[["posterior"]][["forecast"]][["loglik"]]), 2L)
})


test_that("a test sample that does not reach the forecast leaves the object alone", {
  object <- with_forecast(fx_var_fitted())
  scored <- add_predictive_loglik(object, test_sample = object[["data"]][["train"]][["y"]])

  expect_null(scored[["posterior"]][["forecast"]][["loglik"]])
  expect_null(scored[["data"]][["test"]][["y"]])
})


test_that("the log predictive likelihood reaches selection_criteria", {
  object <- with_forecast(fx_var_fitted())
  scored <- add_predictive_loglik(object, test_sample = var_full_data())
  criteria <- selection_criteria(scored)

  expect_true("LPL" %in% names(criteria))
  expect_identical(nrow(attr(criteria[["LPL"]], "terms")), 4L)

  # The criterion is the sum over the periods of the log of the mean of the
  # draws of each, which is what makes the columns add up to a joint density.
  draws <- unclass(scored[["posterior"]][["forecast"]][["loglik"]])
  by_period <- apply(draws, 2, function(x) log(mean(exp(x - max(x)))) + max(x))
  expect_equal(criteria[["LPL"]][["mean"]], sum(by_period))
})


test_that("a model whose states drift is scored under them", {
  object <- with_forecast(fx_var_tvp_fitted("gamma"))
  scored <- add_predictive_loglik(object, test_sample = var_full_data())
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_identical(dim(loglik), c(fx_iterations, 4L))
  expect_true(all(is.finite(loglik)))

  # Holding the states is a different model from letting them drift, so it is a
  # different score. The steps are drawn, hence the seed on either side.
  held <- object
  held[["model"]][["forecast_states"]] <- "hold"
  set.seed(5)
  drifting <- add_predictive_loglik(object, test_sample = var_full_data())
  set.seed(5)
  holding <- add_predictive_loglik(held, test_sample = var_full_data())

  expect_false(isTRUE(all.equal(
    unclass(drifting[["posterior"]][["forecast"]][["loglik"]]),
    unclass(holding[["posterior"]][["forecast"]][["loglik"]]))))
})


test_that("a model whose volatility drifts is scored under it", {
  object <- with_forecast(fx_var_tvp_fitted("sv"))
  scored <- add_predictive_loglik(object, test_sample = var_full_data())
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_identical(dim(loglik), c(fx_iterations, 4L))
  expect_true(all(is.finite(loglik)))
})


test_that("a VEC model is scored against the levels", {
  object <- with_forecast(fx_vec_fitted())
  scored <- add_predictive_loglik(object, test_sample = vec_full_data())
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_s3_class(loglik, "mcmc")
  expect_identical(dim(loglik), c(fx_iterations, 4L))
  expect_true(all(is.finite(loglik)))

  # The levels themselves, not the differences the model was estimated on: the
  # forecast is of the levels and so is what it is scored against.
  starts_at <- stats::tsp(object[["data"]][["train"]][["y"]])[2] + 1 / 4
  expect_equal(as.numeric(scored[["data"]][["test"]][["y"]]),
               as.numeric(stats::window(vec_full_data(), start = starts_at,
                                        end = starts_at + 3 / 4)))

  expect_true("LPL" %in% names(selection_criteria(scored)))
})


test_that("a VEC model whose space drifts is scored under it", {
  object <- with_forecast(fx_vec_tvp_fitted("gamma"))
  scored <- add_predictive_loglik(object, test_sample = vec_full_data())
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]

  expect_identical(dim(loglik), c(fx_iterations, 4L))
  expect_true(all(is.finite(loglik)))
})


test_that("a model that cannot be scored is refused", {
  object <- with_forecast(fx_var_fitted())

  unforecast <- object
  unforecast[["posterior"]][["forecast"]] <- NULL
  expect_error(add_predictive_loglik(unforecast, test_sample = var_full_data()),
               "does not contain forecasts")

  expect_error(add_predictive_loglik(object), "carries none in data\\$test\\$y")

  # The whitelist of algorithms with a normal likelihood is the one the
  # expanding window density uses, and it keeps quantile models out of both.
  quantile_model <- object
  quantile_model[["model"]][["algorithm"]] <- "VarTvpAld"
  expect_error(add_predictive_loglik(quantile_model, test_sample = var_full_data()),
               "not available for the algorithm 'VarTvpAld'")

  wide <- object
  wide[["data"]][["test"]][["y"]] <- cbind(var_full_data(), 1)[1:4, , drop = FALSE]
  expect_error(add_predictive_loglik(wide), "columns, but the model has")

  long <- object
  long[["data"]][["test"]][["y"]] <- as.matrix(var_full_data())
  expect_error(add_predictive_loglik(long), "more than the 4")
})
