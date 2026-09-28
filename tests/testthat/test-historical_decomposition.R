# historical_decomposition(), fevd() with bands, and an i.i.d. variable
# estimated beside a panel not observed whole.

hd_model <- function(iterations = fx_iterations, error = "wishart", tvp = FALSE) {
  model <- create_bvarmodel(var_data(), p = 2, error = error, tvp = tvp,
                            iterations = iterations, burnin = fx_burnin)
  coef <- list(v_i = 1, v_i_det = 0.1)
  if (tvp) {
    coef <- c(coef, list(shape = 3, rate = 1e-4))
  }
  sigma <- switch(error, wishart = list(df = "k", scale = 1),
                  sv = list(mu = 0, v_i = 0.01, shape = 3, rate = 0.01,
                            state_variance = 0.05, offset = 1e-4))
  model <- add_priors(model, coef = coef, sigma = sigma)
  set.seed(31)
  add_posterior_coefficients(add_initial_values(model))
}

test_that("the mean contributions and the baseline add up to the data", {
  model <- hd_model()
  hd <- historical_decomposition(model, response = "Dp")
  expect_s3_class(hd, "bvarhd")
  expect_equal(colnames(hd), c(model$model$endogen, "baseline"))
  expect_equal(nrow(hd), nrow(model$data$train$y))
  expect_equal(stats::tsp(hd), stats::tsp(model$data$train$y))
  expect_equal(as.numeric(rowSums(hd)), as.numeric(model$data$train$y[, "Dp"]), tolerance = 1e-10)
})

test_that("the baseline is the model run from its presample without shocks", {
  model <- hd_model(iterations = 1)
  hd <- historical_decomposition(model, response = "y")
  k <- 3
  p <- 2
  a <- matrix(model$posterior$a$coeffs[1, ], k)
  x <- as.matrix(model$data$train$x)
  tt <- nrow(x)
  # The lags before the sample are the first row's lag columns; every later
  # lag is the baseline's own past.
  b <- matrix(0, tt, k)
  for (t in seq_len(tt)) {
    lags <- x[t, 1:(k * p)]
    for (l in seq_len(p)) {
      if (t > l) {
        lags[(l - 1) * k + 1:k] <- b[t - l, ]
      }
    }
    b[t, ] <- a %*% c(lags, x[t, -(1:(k * p))])
  }
  expect_equal(as.numeric(hd[, "baseline"]), b[, 1], tolerance = 1e-8)
})

test_that("a custom impact matrix names the shocks, and ci gives bands", {
  model <- hd_model()
  impact <- diag(3)
  colnames(impact) <- c("s1", "s2", "s3")
  hd <- historical_decomposition(model, response = "r", type = "custom", impact = impact,
                                 ci = 0.68)
  expect_equal(colnames(hd), c("s1", "s2", "s3", "baseline"))
  expect_true(all(attr(hd, "lower") <= attr(hd, "upper")))
  expect_equal(dim(attr(hd, "lower")), dim(hd))
  median <- historical_decomposition(model, response = "r", statistic = "median")
  expect_equal(dim(median), dim(hd))
})

test_that("what the decomposition cannot describe is refused", {
  expect_error(historical_decomposition(hd_model(), response = "nope"), "endogenous")
  expect_error(historical_decomposition(hd_model(), response = "y", type = "gir"), "type")
  expect_error(historical_decomposition(hd_model(), response = "y", ci = 2), "probability")
  expect_error(historical_decomposition(hd_model(error = "sv"), response = "y"), "every period")
  expect_error(historical_decomposition(hd_model(tvp = TRUE), response = "y"), "constant")
  expect_error(historical_decomposition(hd_model(), response = "y", type = "sign"),
               "add_sign_restrictions")
})

test_that("a completed panel is decomposed draw by draw", {
  data <- var_data()
  data[1:6, "Dp"] <- NA
  model <- create_bvarmodel(data, p = 1, missing = "estimate",
                            iterations = fx_iterations, burnin = fx_burnin)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1), sigma = list(df = "k", scale = 1))
  set.seed(32)
  model <- add_posterior_coefficients(add_initial_values(model))
  hd <- historical_decomposition(model, response = "Dp")
  expect_false(anyNA(hd))
  # The data it decomposes is the mean completed series, observed where it was.
  completed <- colMeans(model$posterior$y$coeffs)[seq(2, length(model$data$train$y), by = 3)]
  expect_equal(as.numeric(attr(hd, "data")), as.numeric(completed), tolerance = 1e-8)
  # An impact function that drops draws still reads the completed panel at the
  # rows of the draws it kept.
  half <- historical_decomposition(model, response = "Dp", type = "custom",
                                   impact = function(draw, i) if (i %% 2 == 0) NULL else diag(3))
  kept <- seq(1, nrow(model$posterior$y$coeffs), by = 2)
  kept_mean <- colMeans(model$posterior$y$coeffs[kept, , drop = FALSE])[seq(2, length(model$data$train$y), by = 3)]
  expect_equal(as.numeric(rowSums(half)), as.numeric(kept_mean), tolerance = 1e-8)

  observed <- !is.na(stats::window(data[, "Dp"], start = stats::start(hd)))
  expect_equal(as.numeric(attr(hd, "data"))[observed],
               as.numeric(stats::window(data[, "Dp"], start = stats::start(hd)))[observed],
               tolerance = 1e-8)
})

test_that("fevd() gives the mean as before, and medians and bands on request", {
  model <- hd_model()
  plain <- fevd(model, response = "y", n_ahead = 4)
  explicit <- fevd(model, response = "y", n_ahead = 4, statistic = "mean")
  expect_identical(unclass(plain), unclass(explicit))
  expect_null(attr(plain, "lower"))

  banded <- fevd(model, response = "y", n_ahead = 4, ci = 0.68)
  expect_equal(unclass(banded)[, ], unclass(plain)[, ])
  lower <- attr(banded, "lower")
  upper <- attr(banded, "upper")
  expect_equal(dim(lower), dim(plain))
  expect_true(all(lower >= -1e-12 & upper <= 1 + 1e-12))
  expect_true(all(lower <= upper))

  median <- fevd(model, response = "y", n_ahead = 4, statistic = "median", ci = 0.68)
  expect_true(all(attr(median, "lower") <= median & median <= attr(median, "upper")))
  expect_error(fevd(model, response = "y", ci = 0.68, max_groups = 2), "cannot be combined")
  expect_error(fevd(model, response = "y", statistic = "mode"), "statistic")
})

test_that("an i.i.d. variable not observed at the start is estimated", {
  # A surprise series, white noise correlated with inflation's errors, that
  # starts twelve quarters into the sample.
  data <- var_data()
  set.seed(33)
  surprise <- stats::ts(stats::rnorm(nrow(data)), start = stats::start(data),
                        frequency = stats::frequency(data))
  data <- cbind(surprise = surprise, data)
  colnames(data) <- c("surprise", colnames(var_data()))
  data[1:12, "surprise"] <- NA
  model <- create_bvarmodel(data, p = 1, iid = "surprise", missing = "estimate",
                            iterations = fx_iterations, burnin = fx_burnin)
  expect_true(all(model$data$train$y[is.na(stats::window(data[, "surprise"],
                                                           start = stats::start(model$data$train$y))),
                                     "surprise"] == 0))
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1), sigma = list(df = "k", scale = 1))
  set.seed(34)
  model <- add_posterior_coefficients(add_initial_values(model))
  k <- 4
  restricted <- seq(1, ncol(model$posterior$a$coeffs), by = k)
  expect_true(all(model$posterior$a$coeffs[, restricted] == 0))
  first <- model$posterior$y$coeffs[, 1]
  expect_gt(stats::sd(first), 0)
  expect_false(anyNA(add_posterior_loglik(model)$posterior$loglik))
})
