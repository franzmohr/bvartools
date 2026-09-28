# The structural quantile VAR: create_bvarmodel(quantile_grid = TRUE).
#
# The numerics are the vendored core's and are tested upstream. What is tested
# here is the plumbing, where a mistake is silent: a level whose block is read
# from the wrong columns, or a quantile that never reaches the sampler, gives a
# model that runs and describes something else.

grid_levels <- c(0.1, 0.25, 0.5, 0.75, 0.9)

grid_model <- function(structural = TRUE) {
  create_bvarmodel(var_data(), p = 1, deterministic = "const",
                   structural = structural, error = "ald",
                   quantile = rev(grid_levels), quantile_grid = TRUE,
                   iterations = fx_iterations, burnin = fx_burnin)
}

grid_prepared <- function(object = grid_model()) {
  add_initial_values(add_priors(object, coef = list(v_i = 1),
                                sigma = list(shape = 3, rate = 0.01)))
}

grid_fitted <- function() {
  cached_fixture("quantile_grid_fitted", {
    set.seed(123456)
    add_posterior_coefficients(grid_prepared())
  })
}

test_that("a grid of quantiles is one model carrying every level", {
  object <- grid_model()

  expect_s3_class(object, "bvarmodel")
  expect_identical(object[["model"]][["algorithm"]], "VarNormalAld")
  # Sorted, which the core requires, and no single quantile beside them.
  expect_identical(object[["model"]][["quantiles"]], grid_levels)
  expect_null(object[["model"]][["quantile"]])
})

test_that("what a grid needs is checked where it is specified", {
  make <- function(...) {
    create_bvarmodel(var_data(), p = 1, deterministic = "const", quantile_grid = TRUE,
                     iterations = fx_iterations, burnin = fx_burnin, ...)
  }
  expect_error(make(error = "ald", quantile = 0.5), "at least two distinct")
  expect_error(make(error = "ald", quantile = c(0.5, 0.5)), "at least two distinct")
  expect_error(make(error = "gamma", quantile = c(0.1, 0.9)), "error = \"ald\"")
  expect_error(make(error = "ald", quantile = c(0.1, 0.9), tvp = TRUE), "constant coefficients")
  expect_error(make(error = "ald", quantile = c(0.1, 0.9), quantile_grid = NA), "quantile_grid")
})

test_that("the draws are stacked level by level", {
  object <- grid_fitted()
  k <- object[["model"]][["k"]]
  nparams <- ncol(object[["data"]][["train"]][["z"]])

  expect_s3_class(object[["posterior"]][["a"]][["coeffs"]], "mcmc")
  expect_identical(dim(object[["posterior"]][["a"]][["coeffs"]]),
                   c(fx_iterations, as.integer(length(grid_levels) * nparams)))
  expect_identical(dim(object[["posterior"]][["u_scale"]][["coeffs"]]),
                   c(fx_iterations, as.integer(length(grid_levels) * k)))
  # The latent precisions are a nuisance of each level's chain and not kept.
  expect_null(object[["posterior"]][["u_sigma_inv"]])
  expect_null(object[["posterior"]][["u_omega_inv"]])
})

test_that("the first level is the chain its single-quantile model draws", {
  single <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                             structural = TRUE, error = "ald", quantile = grid_levels[1],
                             iterations = fx_iterations, burnin = fx_burnin)
  set.seed(123456)
  single <- add_posterior_coefficients(grid_prepared(single))

  first <- split_quantile_grid(grid_fitted())[[1]]
  expect_equal(unclass(first[["posterior"]][["a"]][["coeffs"]]),
               unclass(single[["posterior"]][["a"]][["coeffs"]]), ignore_attr = TRUE)
  expect_equal(unclass(first[["posterior"]][["u_scale"]][["coeffs"]]),
               unclass(single[["posterior"]][["u_scale"]][["coeffs"]]), ignore_attr = TRUE)
})

test_that("each level describes the quantile it was asked for", {
  object <- grid_fitted()
  y <- matrix(t(object[["data"]][["train"]][["y"]]))
  z <- object[["data"]][["train"]][["z"]]
  levels <- split_quantile_grid(object)

  expect_s3_class(levels, "modellist")
  expect_identical(names(levels), paste0("q", grid_levels))
  share_below <- vapply(levels, function(level) {
    mean(y - z %*% colMeans(level[["posterior"]][["a"]][["coeffs"]]) < 0)
  }, numeric(1))
  expect_true(all(abs(share_below - grid_levels) < 0.1))
  # Increasing in the level, which a block read from the wrong columns breaks.
  expect_false(is.unsorted(share_below))
})

test_that("the single levels work where the grid is refused", {
  object <- grid_fitted()
  expect_error(summary(object), "split_quantile_grid")
  expect_error(plot(object), "split_quantile_grid")
  expect_error(irf(object, impulse = "y", response = "r"), "split_quantile_grid")

  level <- split_quantile_grid(object)[[2]]
  expect_identical(level[["model"]][["quantile"]], 0.25)
  expect_output(print(summary(level)), "q = 0.25")
  expect_true(all(is.finite(add_posterior_loglik(level)[["posterior"]][["loglik"]])))
  expect_error(split_quantile_grid(level), "does not hold a grid")
})

test_that("the log likelihood is the density the grid describes", {
  object <- add_posterior_loglik(grid_fitted())
  tt <- nrow(object[["data"]][["train"]][["y"]])

  expect_s3_class(object[["posterior"]][["loglik"]], "mcmc")
  expect_identical(dim(object[["posterior"]][["loglik"]]), c(fx_iterations, tt))
  expect_true(all(is.finite(object[["posterior"]][["loglik"]])))
})

test_that("a structural grid forecasts by simulation", {
  object <- add_forecast_input(grid_fitted(), n_ahead = 4)
  set.seed(1)
  object <- add_posterior_forecasts(object)
  fc <- object[["posterior"]][["forecast"]][["forecasts"]]

  expect_s3_class(fc, "mcmc")
  expect_identical(dim(fc), c(fx_iterations, 12L))
  expect_true(all(is.finite(fc)))
  expect_s3_class(predict(object), "bvarprd")

  # A single quantile still does not.
  single <- split_quantile_grid(grid_fitted())[[3]]
  expect_error(add_posterior_forecasts(add_forecast_input(single, n_ahead = 4)),
               "does not forecast")
})

test_that("a fixed level gives quantile paths and a scenario pins them", {
  object <- add_forecast_input(grid_fitted(), n_ahead = 4)

  # At a fixed level nothing is drawn but the level, so the paths repeat.
  set.seed(1)
  a <- add_posterior_forecasts(object, forecast_quantile = 0.5)
  set.seed(2)
  b <- add_posterior_forecasts(a)
  expect_identical(a[["model"]][["forecast_quantile"]], 0.5)
  expect_equal(unclass(a[["posterior"]][["forecast"]][["forecasts"]]),
               unclass(b[["posterior"]][["forecast"]][["forecasts"]]))

  # The first variable raised by one in the first period: pinned there, and
  # the variables ordered after it respond through the contemporaneous
  # coefficients.
  scenario <- data.frame(period = 1, variable = "y",
                         value = mean(a[["posterior"]][["forecast"]][["forecasts"]][, 1]) + 1)
  s <- add_posterior_forecasts(a, scenario = scenario)
  fa <- unclass(a[["posterior"]][["forecast"]][["forecasts"]])
  fs <- unclass(s[["posterior"]][["forecast"]][["forecasts"]])
  expect_true(all(fs[, 1] == scenario[["value"]]))
  expect_false(isTRUE(all.equal(fa[, 2:3], fs[, 2:3])))
  expect_identical(s[["data"]][["forecast"]][["constraints"]][["variable"]], 1)

  # Zero removes the stored level.
  expect_null(add_posterior_forecasts(s, forecast_quantile = 0)[["model"]][["forecast_quantile"]])

  expect_error(add_posterior_forecasts(object, forecast_quantile = 1), "in \\(0, 1\\)")
  expect_error(add_posterior_forecasts(object, scenario = data.frame(period = 9, variable = 1,
                                                                     value = 0)),
               "from 1 to 4")
  expect_error(add_posterior_forecasts(object, scenario = data.frame(period = 1, variable = "x",
                                                                     value = 0)),
               "does not have")
})

test_that("a grid that is not recursive does not forecast", {
  object <- grid_prepared(grid_model(structural = FALSE))
  set.seed(1)
  object <- add_forecast_input(add_posterior_coefficients(object), n_ahead = 2)
  expect_error(add_posterior_forecasts(object), "structural")
})

test_that("a grid survives a write and read round trip", {
  skip_if_not_installed("hdf5r")
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path))
  object <- add_forecast_input(grid_fitted(), n_ahead = 4)
  object[["model"]][["forecast_quantile"]] <- 0.5
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(as.numeric(restored[["model"]][["quantiles"]]), grid_levels)
  expect_equal(restored[["model"]][["forecast_quantile"]], 0.5)
  expect_equal(unclass(restored[["posterior"]][["a"]][["coeffs"]]),
               unclass(object[["posterior"]][["a"]][["coeffs"]]), ignore_attr = TRUE)
  expect_no_error(add_posterior_forecasts(restored))
})

test_that("a scenario conditions the forecast of a Gaussian VAR too", {
  object <- add_forecast_input(fx_var_fitted(), n_ahead = 3)
  object <- add_posterior_forecasts(object, scenario = data.frame(period = c(1, 2),
                                                                  variable = c(2, 2),
                                                                  value = c(0.5, 0.25)))
  fc <- unclass(object[["posterior"]][["forecast"]][["forecasts"]])
  # Column 2 is the second variable in the first period, column 5 in the second.
  expect_equal(unname(fc[, 2]), rep(0.5, nrow(fc)), tolerance = 1e-8)
  expect_equal(unname(fc[, 5]), rep(0.25, nrow(fc)), tolerance = 1e-8)
})
