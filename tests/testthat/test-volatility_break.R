# The volatility break of VarNormalWishart and VecNormalWishart:
# create_bvarmodel(volatility_break = ) and create_bvecmodel(volatility_break = )
# and sigma$volatility in add_priors(). The sampler is BayesTS's and is tested
# there against its conditionals; what is checked here is that a date reaches it
# as the right period, that its draws come back where they are documented, that
# a break of known size is found, and that window(), vec_to_var() and the model
# file carry it.

# A bivariate VAR(1) whose shocks are six, four and three times their size in
# the three months from August 2006 on, the excess decaying at 0.6 a month after
# them: the break of Lenza and Primiceri (2022) with known scales.
volatility_break_data <- function() {
  set.seed(1)
  tt <- 120
  brk <- 80
  s <- rep(1, tt)
  s[brk:(brk + 2)] <- c(6, 4, 3)
  s[(brk + 3):tt] <- 1 + (3 - 1) * 0.6^(1:(tt - brk - 2))
  a <- matrix(c(0.5, 0.1, 0, 0.4), 2)
  y <- matrix(0, tt, 2)
  for (t in 2:tt) {
    y[t, ] <- a %*% y[t - 1, ] + s[t] * stats::rnorm(2)
  }
  stats::ts(y, start = c(2000, 1), frequency = 12, names = c("a", "b"))
}

volatility_prior <- list(shape = 1, rho_shape1 = 3.04, rho_shape2 = 1.51)

volatility_break_fit <- function() {
  cached_fixture("volatility_break_fit", {
    model <- create_bvarmodel(volatility_break_data(), p = 1, deterministic = "const",
                              volatility_break = c(2006, 8), iterations = 400, burnin = 200)
    model <- add_priors(model, coef = list(v_i = 0.1, v_i_det = 0.01),
                        sigma = list(df = "k", scale = 1, volatility = volatility_prior))
    add_posterior_coefficients(add_seed(add_initial_values(model), 3))
  })
}

volatility_break_vec <- function() {
  cached_fixture("volatility_break_vec", {
    model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                              volatility_break = c(2001, 1), iterations = 60, burnin = 30)
    model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0.1),
                        coint = list(v_i = 0, p_tau_i = 1),
                        sigma = list(df = "k", scale = 1, volatility = volatility_prior))
    add_posterior_coefficients(add_seed(add_initial_values(model), 5))
  })
}

test_that("the break is the period of the estimation sample the date names", {
  data <- volatility_break_data()
  model <- create_bvarmodel(data, p = 1, volatility_break = c(2006, 8))
  y <- model[["data"]][["train"]][["y"]]
  expect_identical(model[["model"]][["volatility_break"]], 79L)
  expect_equal(stats::time(y)[79], 2006 + 7 / 12)
  expect_identical(model[["model"]][["volatility_scales"]], 3L)
  # A single time names the same period. Lag orders share a sample, which
  # starts after the largest of them, and so share the break.
  expect_identical(create_bvarmodel(data, p = 1, volatility_break = 2006 + 7 / 12)
                   [["model"]][["volatility_break"]], 79L)
  models <- create_bvarmodel(data, p = 1:2, volatility_break = c(2006, 8))
  expect_identical(vapply(models, function(m) m[["model"]][["volatility_break"]], 1L),
                   c(78L, 78L))
  expect_identical(create_bvarmodel(data, p = 1, volatility_break = c(2006, 8),
                                    volatility_scales = 1)[["model"]][["volatility_scales"]], 1L)
  expect_null(create_bvarmodel(data, p = 1)[["model"]][["volatility_break"]])
})

test_that("a break is refused where the sampler does not draw it or the date does not fit", {
  data <- volatility_break_data()
  expect_error(create_bvarmodel(data, p = 1, error = "gamma", volatility_break = c(2006, 8)),
               "wishart")
  expect_error(create_bvarmodel(data, p = 1, tvp = TRUE, volatility_break = c(2006, 8)),
               "constant")
  expect_error(create_bvecmodel(vec_data(), p = 2, r = 1, error = "sv",
                                volatility_break = c(2001, 1)), "wishart")
  expect_error(create_bvarmodel(data, p = 1, volatility_break = c(2020, 1)),
               "not a period of the estimation sample")
  # The first period of the data is a lag, not a period of the sample.
  expect_error(create_bvarmodel(data, p = 1, volatility_break = c(2000, 1)),
               "not a period of the estimation sample")
  expect_error(create_bvarmodel(data, p = 1, volatility_break = "2006-08"), "date")
  expect_error(create_bvarmodel(data, p = 1, volatility_break = c(2006, 8),
                                volatility_scales = 0), "at least one")
})

test_that("the prior of the break has to be stated, and only where there is one", {
  model <- create_bvarmodel(volatility_break_data(), p = 1, volatility_break = c(2006, 8))
  coef <- list(v_i = 1, v_i_det = 0.1)
  expect_error(add_priors(model, coef = coef, sigma = list(df = "k", scale = 1)),
               "sigma[$]volatility|element 'volatility'")
  expect_error(add_priors(model, coef = coef,
                          sigma = list(df = "k", scale = 1, volatility = list(shape = 1))),
               "rho_shape1")
  expect_error(add_priors(model, coef = coef,
                          sigma = list(df = "k", scale = 1,
                                       volatility = list(shape = -1, rho_shape1 = 3,
                                                         rho_shape2 = 1.5))),
               "positive")
  expect_identical(add_priors(model, coef = coef,
                              sigma = list(df = "k", scale = 1, volatility = volatility_prior))
                   [["priors"]][["volatility"]], volatility_prior)

  plain <- create_bvarmodel(volatility_break_data(), p = 1)
  expect_error(add_priors(plain, coef = coef,
                          sigma = list(df = "k", scale = 1, volatility = volatility_prior)),
               "volatility_break")
  vec <- create_bvecmodel(vec_data(), p = 2, r = 1)
  expect_error(add_priors(vec, coef = coef, coint = list(v_i = 0, p_tau_i = 1),
                          sigma = list(df = "k", scale = 1, volatility = volatility_prior)),
               "volatility_break")
})

test_that("the scales find a break of known size", {
  fit <- volatility_break_fit()
  draws <- fit[["posterior"]][["volatility"]][["coeffs"]]
  expect_s3_class(draws, "mcmc")
  expect_identical(colnames(draws), c("s1", "s2", "s3", "rho"))
  expect_identical(dim(draws), c(400L, 4L))

  # The shocks were six, four and three times their size: every free scale is
  # well above one, which is where the prior puts its mass, and the decay rate
  # is inside the unit interval.
  medians <- apply(as.matrix(draws), 2, stats::median)
  expect_true(all(medians[c("s1", "s2")] > 3))
  expect_gt(medians[["s3"]], 1.5)
  expect_true(all(as.matrix(draws)[, 1:3] >= 1))
  expect_true(all(draws[, "rho"] > 0 & draws[, "rho"] < 1))
})

test_that("the likelihood, forecasts, window, thinning and model file carry the break", {
  fit <- add_posterior_loglik(volatility_break_fit())
  expect_true(all(is.finite(as.matrix(fit[["posterior"]][["loglik"]]))))
  forecast <- add_posterior_forecasts(add_forecast_input(fit, n_ahead = 3))
  expect_identical(dim(forecast[["posterior"]][["forecast"]][["forecasts"]]), c(400L, 6L))

  # The break moves with the start of a window, is dropped by one that ends
  # before it, and refuses one that starts after it.
  expect_identical(window(fit, start = c(2003, 1))[["model"]][["volatility_break"]], 44L)
  expect_warning(early <- window(fit, end = c(2005, 1)), "ends before the volatility break")
  expect_null(early[["model"]][["volatility_break"]])
  expect_null(early[["priors"]][["volatility"]])
  expect_null(early[["posterior"]][["volatility"]])
  expect_error(window(fit, start = c(2008, 1)), "starts after the volatility break")

  expect_identical(nrow(thin(fit, 2)[["posterior"]][["volatility"]][["coeffs"]]), 200L)

  path <- temp_h5_file()
  write_to_hdf5(fit, filename = path)
  back <- read_model_from_hdf5(path)
  expect_identical(as.integer(back[["model"]][["volatility_break"]]), 79L)
  expect_equal(back[["priors"]][["volatility"]][["rho_shape1"]], 3.04)
  expect_equal(unclass(as.matrix(back[["posterior"]][["volatility"]][["coeffs"]])),
               unclass(as.matrix(fit[["posterior"]][["volatility"]][["coeffs"]])),
               ignore_attr = TRUE)
})

test_that("a VEC model's break reaches its VAR in levels unchanged", {
  fit <- volatility_break_vec()
  y <- fit[["data"]][["train"]][["y"]]
  expect_equal(stats::time(y)[fit[["model"]][["volatility_break"]]], 2001)
  expect_identical(colnames(fit[["posterior"]][["volatility"]][["coeffs"]]),
                   c("s1", "s2", "s3", "rho"))

  # The log likelihood and the forecasts of the VEC model are those of its
  # VAR representation, which reads the break and its draws from it.
  fit <- add_posterior_loglik(fit)
  var <- vec_to_var(fit)
  expect_identical(var[["model"]][["volatility_break"]], fit[["model"]][["volatility_break"]])
  expect_identical(var[["posterior"]][["volatility"]], fit[["posterior"]][["volatility"]])
  expect_equal(as.matrix(add_posterior_loglik(var)[["posterior"]][["loglik"]]),
               as.matrix(fit[["posterior"]][["loglik"]]), ignore_attr = TRUE)
  vec_forecast <- add_posterior_forecasts(add_forecast_input(fit, n_ahead = 2))
  var_forecast <- add_posterior_forecasts(add_forecast_input(var, n_ahead = 2))
  expect_equal(as.matrix(var_forecast[["posterior"]][["forecast"]][["forecasts"]]),
               as.matrix(vec_forecast[["posterior"]][["forecast"]][["forecasts"]]),
               ignore_attr = TRUE)

  path <- temp_h5_file()
  write_to_hdf5(fit, filename = path)
  back <- read_model_from_hdf5(path)
  expect_equal(unclass(as.matrix(back[["posterior"]][["volatility"]][["coeffs"]])),
               unclass(as.matrix(fit[["posterior"]][["volatility"]][["coeffs"]])),
               ignore_attr = TRUE)
})

test_that("the criteria count the scales and the rate they decay at", {

  fit <- add_posterior_loglik(volatility_break_fit())
  tt <- nrow(fit[["data"]][["train"]][["y"]])
  k <- fit[["model"]][["k"]]

  # The count AIC and BIC differ by, which is what each of them penalises with.
  counted <- function(x) {
    sel <- selection_criteria(x)
    (sel[["BIC"]][["mean"]] - sel[["AIC"]][["mean"]]) /
      (log(nrow(x[["data"]][["train"]][["y"]])) - 2)
  }

  # A bivariate VAR(1) with a constant: 2 * (2 + 1) coefficients and 3 of the
  # error covariance, then the three free scales and the decay rate.
  plain <- k * (k * 1 + 1) + k * (k + 1) / 2
  expect_equal(counted(fit), plain + 4)
  expect_equal(.volatility_break_nparams(fit), 4)

  # The same model with the break taken off it counts the coefficients alone.
  # Its deviance is another model's, the break being part of the likelihood, so
  # only the count is comparable.
  bare <- fit
  bare[["model"]][["volatility_break"]] <- NULL
  bare[["model"]][["volatility_scales"]] <- NULL
  expect_equal(counted(bare), plain)
  expect_equal(.volatility_break_nparams(bare), 0)

  # What the four were added to is the penalty and not the deviance: the
  # deviance is still the one of the likelihood at the posterior mean, which
  # scales the covariance of the periods after the break.
  sel <- selection_criteria(fit)
  point <- .posterior_mean_model(fit)
  expect_equal(sel[["AIC"]][["mean"]] - 2 * (plain + 4),
               -2 * sum(add_posterior_loglik(point)[["posterior"]][["loglik"]]))
  expect_equal(sel[["BIC"]][["mean"]] - sel[["AIC"]][["mean"]],
               (log(tt) - 2) * (plain + 4))

  # A scale the sample does not reach is not a parameter the data paid for. A
  # break in the last period has one scale and no decay; one three periods from
  # the end has its three scales and still no decay.
  late <- fit
  late[["model"]][["volatility_break"]] <- as.integer(tt)
  expect_equal(.volatility_break_nparams(late), 1)
  late[["model"]][["volatility_break"]] <- as.integer(tt - 2)
  expect_equal(.volatility_break_nparams(late), 3)
  late[["model"]][["volatility_break"]] <- as.integer(tt - 3)
  expect_equal(.volatility_break_nparams(late), 4)

  # And the VEC model counts it beside the free elements of Pi.
  vec <- add_posterior_loglik(volatility_break_vec())
  k_vec <- vec[["model"]][["k"]]
  k_ect <- ncol(vec[["data"]][["train"]][["w"]])
  n_x <- ncol(vec[["data"]][["train"]][["x"]])
  r <- vec[["model"]][["rank"]]
  expect_equal(counted(vec),
               r * (k_vec + k_ect - r) + k_vec * n_x + k_vec * (k_vec + 1) / 2 + 4)
})
