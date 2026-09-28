# add_prior_options(): the adaptive priors, the stationarity condition and the
# steady-state prior of the vendored BayesTS core, reached from R.

po_model <- function(error = "wishart", p = 2) {
  model <- create_bvarmodel(var_data(), p = p, error = error,
                            iterations = fx_iterations, burnin = fx_burnin)
  sigma <- switch(error, wishart = list(df = "k", scale = 1),
                  gamma = list(shape = 3, rate = 1e-4),
                  sv = list(mu = 0, v_i = 0.01, shape = 3, rate = 0.01,
                            state_variance = 0.05, offset = 1e-4))
  add_priors(model, coef = list(v_i = 1, v_i_det = 0.1), sigma = sigma)
}

test_that("the default Minnesota groups are own lags, other lags and the rest", {
  model <- add_prior_options(po_model(), shrinkage = "minnesota")
  group <- as.numeric(model$priors$a$shrinkage$group)
  k <- 3
  a <- matrix(group, nrow = k)
  # Lags of y1 and y2 in both lag blocks: own on the diagonal.
  expect_equal(a[, 1:3], diag(1, 3) + 2 * (1 - diag(1, 3)))
  expect_equal(a[, 4:6], a[, 1:3])
  expect_equal(a[, 7], rep(0, 3))
  expect_equal(as.numeric(model$priors$a$shrinkage$shape), c(3, 3))
  expect_equal(model$model$shrinkage, "minnesota")

  horseshoe <- add_prior_options(po_model(), shrinkage = "horseshoe")
  expect_equal(as.numeric(horseshoe$priors$a$shrinkage$group), c(rep(1, 18), rep(0, 3)))
  expect_null(horseshoe$priors$a$shrinkage$shape)

  none <- add_prior_options(model, shrinkage = "none")
  expect_null(none$priors$a$shrinkage)
  expect_null(none$model$shrinkage)
})

test_that("the adaptive priors return their scales", {
  set.seed(5)
  minnesota <- add_posterior_coefficients(add_initial_values(
    add_prior_options(po_model(), shrinkage = "minnesota", stationary = TRUE)))
  expect_s3_class(minnesota$posterior$a$shrinkage, "mcmc")
  expect_equal(dim(minnesota$posterior$a$shrinkage), c(fx_iterations, 2))
  expect_true(all(minnesota$posterior$a$shrinkage > 0))

  set.seed(6)
  horseshoe <- add_posterior_coefficients(add_initial_values(
    add_prior_options(po_model("gamma"), shrinkage = "horseshoe")))
  expect_equal(dim(horseshoe$posterior$a$local), c(fx_iterations, 21))
})

test_that("the normal-gamma prior has one global rate per lag order", {
  model <- add_prior_options(po_model(), shrinkage = "normal_gamma")
  group <- matrix(as.numeric(model$priors$a$shrinkage$group), nrow = 3)
  expect_equal(group[, 1:3], matrix(1, 3, 3))
  expect_equal(group[, 4:6], matrix(2, 3, 3))
  expect_equal(group[, 7], rep(0, 3))
  expect_equal(as.numeric(model$priors$a$shrinkage$theta), c(0.1, 0.1))
  expect_equal(as.numeric(model$priors$a$shrinkage$shape), c(0.01, 0.01))

  set.seed(10)
  model <- add_posterior_coefficients(add_initial_values(model))
  expect_equal(dim(model$posterior$a$shrinkage), c(fx_iterations, 2))
  expect_equal(dim(model$posterior$a$local), c(fx_iterations, 21))
  expect_true(all(model$posterior$a$local > 0))

  # Written to a file, theta comes back.
  skip_if_not_installed("hdf5r")
  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(model, file)
  back <- read_model_from_hdf5(file)
  expect_equal(as.numeric(back$priors$a$shrinkage$theta), c(0.1, 0.1))
  expect_equal(back$model$shrinkage, "normal_gamma")
})

test_that("theta of the normal-gamma prior can be drawn", {
  expect_error(add_prior_options(po_model(), shrinkage = list(type = "normal_gamma", theta = 0.1,
                                                              theta_rate = 1)),
               "not both")
  model <- add_prior_options(po_model(), shrinkage = list(type = "normal_gamma", theta_rate = 1))
  expect_null(model[["priors"]][["a"]][["shrinkage"]][["theta"]])
  expect_equal(as.numeric(model$priors$a$shrinkage$theta_rate), c(1, 1))

  set.seed(11)
  model <- add_posterior_coefficients(add_initial_values(model))
  expect_s3_class(model$posterior$a$theta, "mcmc")
  expect_equal(dim(model$posterior$a$theta), c(fx_iterations, 2))
  expect_true(all(model$posterior$a$theta > 0))

  skip_if_not_installed("hdf5r")
  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(model, file)
  back <- read_model_from_hdf5(file)
  expect_equal(as.numeric(back$priors$a$shrinkage$theta_rate), c(1, 1))
  expect_equal(unclass(back$posterior$a$theta), unclass(model$posterior$a$theta),
               ignore_attr = TRUE)
})

test_that("every kept draw is stationary under the stationarity condition", {
  set.seed(7)
  model <- add_posterior_coefficients(add_initial_values(
    add_prior_options(po_model(p = 1), stationary = TRUE)))
  a <- model$posterior$a$coeffs
  roots <- apply(a, 1, function(d) max(Mod(eigen(matrix(d, 3)[, 1:3])$values)))
  expect_true(all(roots < 1))
})

test_that("under a steady-state prior the intercept is (I - A) mu in every draw", {
  set.seed(8)
  model <- add_posterior_coefficients(add_initial_values(
    add_prior_options(po_model("sv", p = 1), steady_state = list(mu = c(0.5, 0.5, 0), v_i = 1))))
  expect_equal(dim(model$posterior$mu$coeffs), c(fx_iterations, 3))
  for (d in c(1, fx_iterations)) {
    a <- matrix(model$posterior$a$coeffs[d, ], 3)
    mu <- as.numeric(model$posterior$mu$coeffs[d, ])
    expect_equal(a[, 4], as.numeric((diag(3) - a[, 1:3]) %*% mu), tolerance = 1e-10)
  }
})

test_that("options a model cannot honour are refused", {
  expect_error(add_prior_options(po_model(), shrinkage = "lasso"), "'horseshoe' or 'normal_gamma'")
  expect_error(add_posterior_coefficients(add_initial_values(
    add_prior_options(po_model(), shrinkage = list(type = "normal_gamma", theta = -1)))),
    "theta must be positive")
  expect_error(add_prior_options(po_model(), shrinkage = list(type = "minnesota", group = 1)),
               "one element per coefficient")
  expect_error(add_prior_options(po_model(), steady_state = list(mu = 0)), "'mu' and 'v_i'")
  expect_error(add_prior_options(po_model(), constraints = list(shape = 1, rate = 1)),
               "no soft constraints")
  expect_error(add_prior_options(create_bvarmodel(var_data())), "add_priors")

  # Refused by the core, naming the algorithm.
  tvp <- create_bvarmodel(var_data(), p = 1, tvp = TRUE, iterations = 10, burnin = 0)
  tvp <- add_priors(tvp, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                    sigma = list(df = "k", scale = 1))
  tvp <- add_initial_values(add_prior_options(tvp, stationary = TRUE))
  expect_error(add_posterior_coefficients(tvp), "VarTvpWishart")
})

test_that("the options survive the round trip through a file", {
  skip_if_not_installed("hdf5r")
  model <- add_initial_values(add_prior_options(po_model(), shrinkage = "minnesota", stationary = TRUE))
  set.seed(9)
  model <- add_posterior_coefficients(model)
  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(model, file)
  back <- read_model_from_hdf5(file)
  expect_equal(back$model$shrinkage, "minnesota")
  expect_true(back$model$stationary)
  expect_equal(back$priors$a$shrinkage[names(model$priors$a$shrinkage)],
               model$priors$a$shrinkage, ignore_attr = TRUE)
  expect_equal(unclass(back$posterior$a$shrinkage), unclass(model$posterior$a$shrinkage),
               ignore_attr = TRUE)
})
