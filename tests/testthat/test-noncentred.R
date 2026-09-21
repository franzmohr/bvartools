# The non-centred prior 'omega_v' of a VAR with time varying parameters and
# stochastic volatility, and the draws of the Savage-Dickey test for time
# variation it brings (Chan 2018). The sampler itself is the vendored BayesTS
# core, which tests the numbers upstream -- that the ordinates are the
# densities they claim to be and that the Bayes factors come out on the side
# simulated data were generated on. What is tested here is the R side: that
# add_priors() writes the prior, refuses it where no sampler reads it, and that
# the draws come back, survive a file and leave the log likelihood computable.

nc_sigma_prior <- list(mu = 0, v_i = 0.01, omega_v = 0.1, state_variance = 0.05,
                       offset = 1e-8)

nc_model <- function(error = "sv+covar", coef = list(v_i = 1, v_i_det = 0.1, omega_v = 0.001),
                     sigma = nc_sigma_prior) {
  model <- create_bvarmodel(at_data_var(), p = 1, deterministic = "const", tvp = TRUE,
                            error = error, iterations = fx_iterations, burnin = fx_burnin)
  add_priors(model, coef = coef, sigma = sigma)
}

# Output growth, inflation and the change in the interest rate of at_macrodata,
# in percent, over the whole sample.
at_data_var <- function() {
  levels <- at_domestic()
  stats::ts.intersect(y = diff(levels[, "y"]), Dp = levels[, "Dp"],
                      r = diff(levels[, "r"])) * 100
}

test_that("add_priors() writes omega_v in place of shape and rate", {
  model <- nc_model()
  k <- model[["model"]][["k"]]
  n_a <- nrow(model[["priors"]][["a"]][["mu"]])

  expect_equal(as.numeric(model[["priors"]][["a"]][["omega_v"]]), rep(0.001, n_a))
  expect_equal(as.numeric(model[["priors"]][["psi"]][["omega_v"]]), rep(0.001, k * (k - 1) / 2))
  expect_equal(as.numeric(model[["priors"]][["u_sigma"]][["omega_v"]]), rep(0.1, k))
  for (block in c("a", "psi", "u_sigma")) {
    expect_null(model[["priors"]][[block]][["shape"]])
    expect_null(model[["priors"]][[block]][["rate"]])
  }

  # The state variances start at their prior mean, omega_v.
  model <- add_initial_values(model)
  expect_equal(diag(model[["initial"]][["a_sigma_inv"]]), rep(1 / 0.001, n_a))
})

test_that("each block chooses its prior for itself", {
  model <- nc_model(coef = list(v_i = 1, shape = 3, rate = 0.01))
  expect_null(model[["priors"]][["a"]][["omega_v"]])
  expect_false(is.null(model[["priors"]][["a"]][["shape"]]))
  expect_equal(as.numeric(model[["priors"]][["u_sigma"]][["omega_v"]]), rep(0.1, 3))
})

test_that("omega_v is refused where no sampler reads it, and beside shape and rate", {
  expect_error(nc_model(coef = list(v_i = 1, omega_v = 0.001, shape = 3, rate = 0.01)),
               "both 'omega_v'")
  expect_error(nc_model(sigma = c(nc_sigma_prior, shape = 3, rate = 0.01)), "both 'omega_v'")
  expect_error(nc_model(coef = list(v_i = 1, omega_v = -1)), "finite positive")
  expect_error(nc_model(coef = list(v_i = 1, omega_v = c(0.1, 0.2))), "finite positive")

  gamma <- create_bvarmodel(at_data_var(), p = 1, deterministic = "const", tvp = TRUE,
                            error = "gamma", iterations = 10, burnin = 10)
  expect_error(add_priors(gamma, coef = list(v_i = 1, omega_v = 0.001),
                          sigma = list(shape = 3, rate = 0.01)),
               "only available")

  constant <- create_bvarmodel(at_data_var(), p = 1, deterministic = "const", tvp = FALSE,
                               error = "sv", iterations = 10, burnin = 10)
  expect_error(add_priors(constant, coef = list(v_i = 1), sigma = nc_sigma_prior),
               "only available")
})

test_that("the draws of the test come back beside sigma and survive a file", {
  model <- add_initial_values(nc_model())
  model <- add_posterior_coefficients(add_seed(model, 211))
  k <- model[["model"]][["k"]]

  for (block in c("a", "psi", "u_sigma_inv")) {
    draws <- model[["posterior"]][[block]]
    expect_true(all(c("omega", "omega_log_zero", "omega_log_zero_joint") %in% names(draws)))
    expect_s3_class(draws[["omega"]], "mcmc")
    expect_equal(ncol(draws[["omega"]]), ncol(draws[["sigma"]]))
    expect_equal(ncol(draws[["omega_log_zero_joint"]]), 1)
    # sigma is still the state variance, now the square of omega.
    expect_equal(as.numeric(draws[["sigma"]]), as.numeric(draws[["omega"]])^2)
    expect_true(all(is.finite(draws[["omega_log_zero"]])))
  }
  # The log-volatilities are independent given the draw under a diagonal
  # v_inv, so their joint ordinate is the sum of the per-state ones.
  expect_equal(as.numeric(model[["posterior"]][["u_sigma_inv"]][["omega_log_zero_joint"]]),
               rowSums(as.matrix(model[["posterior"]][["u_sigma_inv"]][["omega_log_zero"]])))

  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(model, file)
  back <- read_model_from_hdf5(file)
  expect_equal(as.numeric(back[["priors"]][["u_sigma"]][["omega_v"]]), rep(0.1, k))
  expect_equal(as.numeric(back[["priors"]][["a"]][["omega_v"]]),
               as.numeric(model[["priors"]][["a"]][["omega_v"]]))
  expect_equal(as.matrix(back[["posterior"]][["a"]][["omega_log_zero"]]),
               as.matrix(model[["posterior"]][["a"]][["omega_log_zero"]]),
               ignore_attr = TRUE)
})

test_that("a non-centred posterior is scored like any other", {
  model <- add_initial_values(nc_model())
  model <- add_posterior_coefficients(add_seed(model, 212))
  model <- add_posterior_loglik(model)
  expect_true(all(is.finite(model[["posterior"]][["loglik"]])))
})

test_that("time_variation_test() reports every state and block of the prior", {
  model <- add_initial_values(nc_model())
  model <- add_posterior_coefficients(add_seed(model, 213))
  k <- model[["model"]][["k"]]
  n_a <- nrow(model[["priors"]][["a"]][["mu"]])
  n_psi <- k * (k - 1) / 2

  res <- time_variation_test(model)
  expect_s3_class(res, "bvartimevar")
  expect_equal(nrow(res), n_a + n_psi + k + 3)
  expect_equal(as.vector(table(res[["block"]])[c("coefficients", "covariances", "volatilities")]),
               c(n_a + 1, n_psi + 1, k + 1))
  expect_true(all(is.finite(res[["log_bf"]])))
  expect_true(all(res[["nse"]] >= 0))

  # The coefficients are vec(A): the equation cycles fastest.
  coefs <- res[res[["block"]] == "coefficients" & res[["term"]] != "(joint)", ]
  expect_equal(coefs[["equation"]][seq_len(k)], model[["model"]][["endogen"]])
  expect_equal(unique(coefs[["term"]][seq_len(k)]), coefs[["term"]][1])
  # Psi row by row: (2,1), (3,1), (3,2).
  y <- model[["model"]][["endogen"]]
  covs <- res[res[["block"]] == "covariances" & res[["term"]] != "(joint)", ]
  expect_equal(covs[["equation"]], y[c(2, 3, 3)])
  expect_equal(covs[["term"]], y[c(1, 1, 2)])

  # The Savage-Dickey ratio, spelled out for one log-volatility.
  lz <- as.numeric(model[["posterior"]][["u_sigma_inv"]][["omega_log_zero"]][, 2])
  expected <- dnorm(0, 0, sqrt(0.1), log = TRUE) - log(mean(exp(lz)))
  vol <- res[res[["block"]] == "volatilities" & res[["term"]] != "(joint)", ]
  expect_equal(vol[["log_bf"]][2], expected)

  expect_equal(nrow(time_variation_test(model, joint = FALSE)), n_a + n_psi + k)
  expect_output(print(res), "Volatilities")

  # Thinning keeps the draws of the test with the rest of the posterior.
  thinned <- thin(model, thin = 2)
  expect_equal(nrow(thinned[["posterior"]][["a"]][["omega_log_zero"]]),
               nrow(thinned[["posterior"]][["a"]][["coeffs"]]))
})

test_that("time_variation_test() reports only the blocks drawn under omega_v", {
  model <- add_initial_values(nc_model(coef = list(v_i = 1, shape = 3, rate = 0.01)))
  model <- add_posterior_coefficients(add_seed(model, 214))
  res <- time_variation_test(model)
  expect_equal(unique(res[["block"]]), "volatilities")

  centred <- add_initial_values(nc_model(coef = list(v_i = 1, shape = 3, rate = 0.01),
                                         sigma = tvp_sigma_prior("sv")))
  centred <- add_posterior_coefficients(add_seed(centred, 215))
  expect_error(time_variation_test(centred), "non-centred prior")
  expect_error(time_variation_test(list()), "bvarmodel")
})
