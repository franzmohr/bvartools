# Quantile VARs: create_bvarmodel(error = "ald") through to posterior draws.
#
# The failure mode these models have and the others do not is silent. A model
# that has lost its quantile somewhere estimates the median instead and looks
# entirely healthy, so the tests below check the estimand itself -- the share of
# fitted residuals below zero -- rather than only the shape of the output.

ald_model <- function(quantile = 0.25, tvp = FALSE) {
  create_bvarmodel(var_data(), p = 1, deterministic = "const",
                   error = "ald", quantile = quantile, tvp = tvp,
                   iterations = fx_iterations, burnin = fx_burnin)
}

ald_fitted <- function(quantile = 0.25) {
  cached_fixture(paste0("ald_fitted_", quantile), {
    object <- add_priors(ald_model(quantile),
                         coef = list(v_i = 1),
                         sigma = list(shape = 3, rate = 0.01))
    set.seed(987654)
    add_posterior_coefficients(add_initial_values(object))
  })
}

test_that("a quantile is part of the specification", {
  object <- ald_model(quantile = 0.8)

  expect_s3_class(object, "bvarmodel")
  expect_identical(object[["model"]][["algorithm"]], "VarNormalAld")
  expect_identical(object[["model"]][["error"]], "ald")
  expect_identical(object[["model"]][["quantile"]], 0.8)
  expect_identical(ald_model(tvp = TRUE)[["model"]][["algorithm"]], "VarTvpAld")
})

test_that("a grid of quantiles is a list of models", {
  objects <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                              error = "ald", quantile = c(0.1, 0.5, 0.9),
                              iterations = fx_iterations, burnin = fx_burnin)

  expect_s3_class(objects, "modellist")
  expect_length(objects, 3)
  expect_identical(vapply(objects, function(x) {x[["model"]][["quantile"]]}, numeric(1)),
                   c(0.1, 0.5, 0.9))
  # Only the quantile differs; the data of the three is the same.
  expect_identical(objects[[1]][["data"]], objects[[3]][["data"]])
})

test_that("the quantile has to be a probability", {
  expect_error(ald_model(quantile = 0), "between 0 and 1")
  expect_error(ald_model(quantile = 1.5), "between 0 and 1")
  expect_error(ald_model(quantile = NA_real_), "non-missing")
  # Every other model ignores the argument and does not carry the field.
  expect_null(create_bvarmodel(var_data(), p = 1, quantile = 0.2,
                               iterations = fx_iterations,
                               burnin = fx_burnin)[["model"]][["quantile"]])
})

test_that("the prior and the initial values cover the asymmetric Laplace", {
  object <- add_priors(ald_model(), coef = list(v_i = 1),
                       sigma = list(shape = 3, rate = 0.01))
  k <- object[["model"]][["k"]]

  # A scale per equation, not an error variance-covariance matrix.
  expect_identical(dim(object[["priors"]][["u_scale"]][["shape"]]), c(k, 1L))
  expect_identical(dim(object[["priors"]][["u_scale"]][["rate"]]), c(k, 1L))
  expect_null(object[["priors"]][["u_sigma"]])

  object <- add_initial_values(object)
  tt <- nrow(object[["data"]][["train"]][["y"]])
  expect_identical(dim(object[["initial"]][["w"]]), c(tt, k))
  expect_true(all(object[["initial"]][["w"]] > 0))
  expect_identical(dim(object[["initial"]][["u_scale"]]), c(k, 1L))
  expect_true(all(object[["initial"]][["u_scale"]] > 0))
  expect_no_error(.check_bvarpost_input(object))
})

test_that("posterior draws arrive in the layout the package uses", {
  object <- ald_fitted()
  k <- object[["model"]][["k"]]
  tt <- nrow(object[["data"]][["train"]][["y"]])

  expect_s3_class(object[["posterior"]][["a"]][["coeffs"]], "mcmc")
  expect_identical(dim(object[["posterior"]][["a"]][["coeffs"]]), c(fx_iterations, 12L))
  # One scale per equation and, like the stochastic volatility models, a
  # precision per period.
  expect_identical(dim(object[["posterior"]][["u_scale"]][["coeffs"]]),
                   c(fx_iterations, k))
  expect_true(all(object[["posterior"]][["u_scale"]][["coeffs"]] > 0))
  expect_identical(dim(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
                   c(fx_iterations, as.integer(k * k * tt)))
})

test_that("the coefficients describe the quantile they were asked for", {
  y <- matrix(ald_fitted(0.25)[["data"]][["train"]][["y"]])
  z <- ald_fitted(0.25)[["data"]][["train"]][["z"]]

  share_below <- function(quantile) {
    a <- colMeans(ald_fitted(quantile)[["posterior"]][["a"]][["coeffs"]])
    mean(y - z %*% a < 0)
  }

  # The defining property of the estimand, and the one thing a model that has
  # dropped its skew term somewhere would fail: the share of residuals below
  # zero is the quantile. Few draws, so the tolerance is wide -- a model
  # estimating the median instead would come out near 0.5 for both.
  expect_equal(share_below(0.25), 0.25, tolerance = 0.1)
  expect_equal(share_below(0.75), 0.75, tolerance = 0.1)
})

test_that("the log likelihood is added like any other model's", {
  object <- add_posterior_loglik(ald_fitted())
  tt <- nrow(object[["data"]][["train"]][["y"]])

  expect_s3_class(object[["posterior"]][["loglik"]], "mcmc")
  expect_identical(dim(object[["posterior"]][["loglik"]]), c(fx_iterations, tt))
  expect_true(all(is.finite(object[["posterior"]][["loglik"]])))
})

test_that("time varying coefficients draw a path per period", {
  object <- add_priors(ald_model(quantile = 0.5, tvp = TRUE),
                       coef = list(v_i = 1, shape = 3, rate = 1e-04),
                       sigma = list(shape = 3, rate = 0.01))
  set.seed(987654)
  object <- add_posterior_coefficients(add_initial_values(object))
  tt <- nrow(object[["data"]][["train"]][["y"]])

  expect_identical(dim(object[["posterior"]][["a"]][["coeffs"]]),
                   c(fx_iterations, as.integer(12 * tt)))
  expect_identical(dim(object[["posterior"]][["a"]][["sigma"]]), c(fx_iterations, 12L))
  expect_no_error(add_posterior_loglik(object))
})

test_that("the summary reports the quantile and the applications still work", {
  object <- add_posterior_loglik(ald_fitted())

  expect_s3_class(summary(object), "summary.bvarmodel")
  expect_output(print(summary(object)), "Quantile-VAR")
  expect_output(print(summary(object)), "q = 0.25")
  expect_s3_class(irf(object, impulse = "income", response = "cons", n_ahead = 3),
                  "bvarirf")
  expect_s3_class(fevd(object, response = "cons", n_ahead = 3), "bvarfevd")
})

test_that("what a quantile model does not do is refused with a reason", {
  # No forecast: the h step quantile is not the quantile of the iterated one
  # step quantiles.
  expect_error(add_posterior_forecasts(ald_fitted()),
               "does not forecast")
  # No SSVS, refused when the specification is made.
  expect_error(create_bvarmodel(var_data(), p = 1, deterministic = "const",
                                error = "ald", varsel = "ssvs",
                                iterations = fx_iterations, burnin = fx_burnin),
               "not available for a quantile")
  # And again when a specification is edited to ask for it after the fact,
  # which is the only way past the refusal above.
  edited <- ald_model()
  edited[["model"]][["varsel"]] <- "ssvs"
  expect_error(add_priors(edited, coef = list(v_i = 1),
                          sigma = list(shape = 3, rate = 0.01),
                          varsel = list(inprior = 0.5, tau = c(0.05, 10))),
               "SSVS is not available")
})
