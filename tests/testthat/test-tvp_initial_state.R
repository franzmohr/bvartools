# The level of a time varying coefficient path comes from the data, not from
# the starting values.
#
# Regression test for the vendored BayesTS core. Its time varying samplers drew
# the path against the previous draw of the state before the sample, with the
# random walk's innovation variance as the prior covariance of the first period,
# and with a variance as small as a rate of 1e-12 makes it the chain returned
# the output of add_initial_values() as the posterior. Two chains started five
# units apart would differ by five; the tolerance below is a fraction of that.

tvp_sigma_prior <- function(error) {
  if (error == "sv") {
    list(mu = 0, v_i = 0.01, shape = 3, rate = 0.01,
         state_variance = 0.05, offset = 1e-8)
  } else {
    list(shape = 3, rate = 0.01)
  }
}

test_that("time varying coefficients do not stay at their starting values", {
  for (error in c("gamma", "sv")) {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", tvp = TRUE,
                              error = error, iterations = 200, burnin = 200)
    model <- add_priors(model,
                        coef = list(v_i = 1, v_i_det = 0.1, shape = 3,
                                    rate = 1e-12, rate_det = 1e-12),
                        sigma = tvp_sigma_prior(error))
    set.seed(1)
    model <- add_initial_values(model)

    far <- model
    far[["initial"]][["a"]][] <- far[["initial"]][["a"]] + 5
    far[["initial"]][["a_init"]][] <- far[["initial"]][["a_init"]] + 5

    set.seed(2)
    near <- add_posterior_coefficients(model)
    set.seed(2)
    far <- add_posterior_coefficients(far)

    # Posterior mean of each coefficient, over the draws and the periods.
    nparams <- ncol(model[["data"]][["train"]][["z"]])
    path_mean <- function(object) {
      rowMeans(matrix(colMeans(object[["posterior"]][["a"]][["coeffs"]]), nparams))
    }

    expect_lt(max(abs(path_mean(near) - path_mean(far))), 1, label = error)
  }
})

test_that("a time varying model needs a proper prior on the state before the sample", {
  # The draw integrates that state out of the prior of the first period, which
  # takes the inverse of its prior precision.
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", tvp = TRUE,
                            error = "gamma", iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0, shape = 3, rate = 1e-4),
                      sigma = list(shape = 3, rate = 0.01))
  set.seed(1)
  model <- suppressWarnings(add_initial_values(model))

  expect_error(add_posterior_coefficients(model), "positive definite")
})
