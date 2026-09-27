# A panel not observed whole: create_bvarmodel(missing = "estimate",
# aggregate, soft), what it writes to data$train$constraints and what the
# samplers make of it. The numerics are the vendored BayesTS core's and tested
# there; what is tested here is the translation.

# var_data() with the output growth of a quarter observed only as the mean of
# it and the two before it, every third quarter, and one inflation rate lost.
mf_data <- function() {
  data <- var_data()
  mean3 <- stats::filter(data[, "y"], rep(1 / 3, 3), sides = 1)
  data[, "y"] <- ifelse(seq_len(nrow(data)) %% 3 == 0, mean3, NA)
  data[10, "Dp"] <- NA
  data
}

mf_priors <- function(model) {
  add_priors(model, coef = list(v_i = 1, v_i_det = 0.1), sigma = list(df = "k", scale = 1))
}

test_that("aggregation_weights() returns the weights it documents", {
  expect_equal(aggregation_weights("average", 3), rep(1 / 3, 3))
  expect_equal(aggregation_weights("sum", 4), rep(1, 4))
  expect_equal(aggregation_weights("growth", 3), c(1, 2, 3, 2, 1) / 3)
  expect_error(aggregation_weights("sum", 0), "whole number")
})

test_that("missing = 'omit' leaves a complete panel as it was", {
  data <- var_data()
  expect_identical(create_bvarmodel(data, p = 1),
                   create_bvarmodel(data, p = 1, missing = "estimate"))
})

test_that("every observation becomes one constraint row, placed by date", {
  data <- mf_data()
  # The sample starts in the second quarter, so the first average, in the
  # third, would reach into the first.
  expect_message(
    model <- create_bvarmodel(data, p = 1, iterations = 10, burnin = 0,
                              aggregate = list(y = aggregation_weights("average", 3))),
    "1 observation(s) of an aggregated series reach before the estimation sample", fixed = TRUE)
  c <- model$data$train$constraints
  y <- model$data$train$y
  expect_named(c, c("value", "group", "row", "period", "variable", "weight"))
  expect_equal(length(c$value), sum(!is.na(stats::window(data, start = stats::tsp(y)[1]))) - 1)
  expect_true(all(c$group == 0))

  # A single observation pins its entry, and y agrees with it there.
  single <- which(tabulate(c$row) == 1)
  entries <- match(single, c$row)
  expect_equal(c$value[single],
               as.numeric(y)[(c$variable[entries] - 1) * nrow(y) + c$period[entries]])

  # An aggregate reaches over the three quarters it averages, oldest first,
  # and its value is that of the quarter it ends in.
  first <- which(c$variable == 1)[1:3]
  expect_equal(c$weight[first], rep(1 / 3, 3))
  expect_equal(diff(c$period[first]), c(1, 1))
  row <- c$row[first[1]]
  end_time <- stats::time(y)[c$period[first[3]]]
  expect_equal(c$value[row], as.numeric(stats::window(data[, "y"], start = end_time, end = end_time)))

  # The gap is not a hole in the regressors.
  expect_false(anyNA(model$data$train$z))
  # And the data the model was made from stay as they were given.
  expect_identical(model$data$original$endogen, data)
})

test_that("arguments that cannot be honoured are refused", {
  data <- mf_data()
  expect_error(create_bvarmodel(data, missing = "drop"), "'omit' or 'estimate'")
  expect_error(create_bvarmodel(data, aggregate = list(nope = 1)), "not in 'data'")
  expect_error(create_bvarmodel(data, aggregate = list(1)), "named list")
  expect_error(create_bvarmodel(data, soft = "y"), "missing = 'estimate'")
  expect_error(create_bvarmodel(data, missing = "estimate", error = "ald"), "not available")
  expect_error(create_bvarmodel(data, missing = "estimate", structural = TRUE, error = "gamma"),
               "not available")
})

test_that("the samplers complete the panel and score what was observed", {
  for (error in c("wishart", "gamma", "sv")) {
    for (tvp in c(FALSE, TRUE)) {
      model <- create_bvarmodel(mf_data(), p = 1, error = error, tvp = tvp,
                                iterations = fx_iterations, burnin = fx_burnin,
                                aggregate = list(y = aggregation_weights("average", 3)))
      sigma <- switch(error, wishart = list(df = "k", scale = 1),
                      gamma = list(shape = 3, rate = 1e-4),
                      sv = list(mu = 0, v_i = 0.01, shape = 3, rate = 0.01,
                                state_variance = 0.05, offset = 1e-4))
      coef <- list(v_i = 1, v_i_det = 0.1)
      if (tvp) {
        coef <- c(coef, list(shape = 3, rate = 1e-4))
      }
      model <- add_priors(model, coef = coef, sigma = sigma)
      model <- add_initial_values(model)
      set.seed(1)
      model <- add_posterior_coefficients(model)
      model <- add_posterior_loglik(model)

      y <- model$data$train$y
      label <- paste(error, tvp)
      expect_s3_class(model$posterior$y$coeffs, "mcmc")
      expect_equal(dim(model$posterior$y$coeffs), c(fx_iterations, length(y)), label = label)
      expect_equal(ncol(model$posterior$loglik), nrow(y), label = label)
      expect_false(anyNA(model$posterior$loglik), label = label)

      # Every draw honours every hard row: y is stacked by period.
      c <- model$data$train$constraints
      k <- ncol(y)
      draw <- as.numeric(model$posterior$y$coeffs[fx_iterations, ])
      fitted <- tapply(c$weight * draw[(c$period - 1) * k + c$variable], c$row, sum)
      expect_equal(as.numeric(fitted), c$value, tolerance = 1e-8, label = label)
    }
  }
})

test_that("a forecast starts from each draw's completed panel", {
  model <- create_bvarmodel(mf_data(), p = 1, iterations = fx_iterations, burnin = fx_burnin,
                            aggregate = list(y = aggregation_weights("average", 3)))
  model <- add_initial_values(mf_priors(model))
  set.seed(2)
  model <- add_posterior_coefficients(model)
  model <- add_forecast_input(model)
  model <- add_posterior_forecasts(model)
  expect_false(anyNA(model$posterior$forecast$forecasts))
})

test_that("soft rows draw their error precision, and need its prior", {
  model <- create_bvarmodel(mf_data(), p = 1, iterations = fx_iterations, burnin = fx_burnin,
                            aggregate = list(y = aggregation_weights("average", 3)), soft = "y")
  c <- model$data$train$constraints
  expect_true(all(c$group[unique(c$row[c$variable == 1])] == 1))
  expect_true(all(c$group[unique(c$row[c$variable != 1])] == 0))

  model <- add_initial_values(mf_priors(model))
  expect_error(add_posterior_coefficients(model), "add_prior_options")

  model <- add_prior_options(model, constraints = list(shape = 2, rate = 0.1))
  expect_equal(as.numeric(model$priors$constraints$shape), 2)
  set.seed(3)
  model <- add_posterior_coefficients(model)
  expect_equal(dim(model$posterior$constraints_inv$coeffs), c(fx_iterations, 1))
  expect_true(all(model$posterior$constraints_inv$coeffs > 0))
})

test_that("a panel not observed whole survives the round trip through a file", {
  skip_if_not_installed("hdf5r")
  model <- create_bvarmodel(mf_data(), p = 1, iterations = fx_iterations, burnin = fx_burnin,
                            aggregate = list(y = aggregation_weights("average", 3)), soft = "y")
  model <- add_prior_options(mf_priors(model))
  model <- add_initial_values(model)
  set.seed(4)
  model <- add_posterior_coefficients(model)

  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(model, file)
  back <- read_model_from_hdf5(file)

  for (i in names(model$data$train$constraints)) {
    expect_equal(back$data$train$constraints[[i]], as.numeric(model$data$train$constraints[[i]]),
                 label = i)
  }
  expect_equal(back$priors$constraints[c("shape", "rate")], model$priors$constraints,
               ignore_attr = TRUE)
  expect_equal(unclass(back$posterior$y$coeffs), unclass(model$posterior$y$coeffs),
               ignore_attr = TRUE)
  expect_equal(unclass(back$posterior$constraints_inv$coeffs),
               unclass(model$posterior$constraints_inv$coeffs), ignore_attr = TRUE)

  # What was read back scores as the original does.
  expect_equal(add_posterior_loglik(back)$posterior$loglik,
               add_posterior_loglik(model)$posterior$loglik, ignore_attr = TRUE)
})
