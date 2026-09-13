# The autocorrelation of the state equation of the cointegration space.
#
# Koop, Leon-Gonzalez and Strachan (2011) treat rho as a parameter rather than a
# fixed hyperparameter. Giving `coint` the support of a uniform prior on it --
# `rho_min` and `rho_max` -- turns the draw on in all three algorithms whose
# cointegration space moves, and `coint$rho` becomes the value the chain starts
# at. Without the pair nothing changes and the posterior carries no rho.

vec_tvp_fitted <- function(error, coint, iterations = 20, burnin = 10) {

  data("us_macrodata", envir = environment())

  object <- create_bvecmodel(data = us_macrodata, p = 2, const = "unrestricted",
                             r = 1, tvp = TRUE, error = error,
                             iterations = iterations, burnin = burnin)

  sigma_prior <- switch(error,
                        sv = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                  state_variance = 0.05, offset = 1e-8),
                        gamma = list(shape = 3, rate = 0.01),
                        list(df = 3, scale = 1))

  object <- add_priors(object,
                       coef = list(v_i = 1, v_i_det = 1 / 10, shape = 3,
                                   rate = 1e-8, rate_det = 1e-8),
                       sigma = sigma_prior,
                       coint = coint)
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

test_that("every time varying cointegration algorithm can draw rho", {
  specifications <- list(
    list(error = "wishart", algorithm = "VecTvpWishart"),
    list(error = "gamma",   algorithm = "VecTvpGamma"),
    list(error = "sv",      algorithm = "VecTvpStochvol")
  )

  for (specification in specifications) {

    set.seed(456789)
    fixed <- vec_tvp_fitted(specification[["error"]], list(rho = 0.999))

    set.seed(456789)
    drawn <- vec_tvp_fitted(specification[["error"]],
                            list(rho = 0.99, rho_min = 0.9, rho_max = 0.999))

    expect_equal(drawn[["model"]][["algorithm"]], specification[["algorithm"]])

    # Held fixed, rho is a hyperparameter the caller already has in the prior,
    # so handing a column of the same number back would read as a posterior.
    expect_null(fixed[["posterior"]][["beta"]][["rho"]])

    rho <- drawn[["posterior"]][["beta"]][["rho"]]
    expect_s3_class(rho, "mcmc")
    expect_equal(dim(rho), c(20L, 1L))

    # Inside the support it was given, and moving within it -- a draw that never
    # moves is what a block that silently did nothing would leave behind.
    expect_true(all(rho >= 0.9))
    expect_true(all(rho <= 0.999))
    expect_gt(length(unique(as.numeric(rho))), 1)

    # The cointegration path is still there and is still the full one.
    expect_equal(ncol(drawn[["posterior"]][["beta"]][["coeffs"]]),
                 ncol(fixed[["posterior"]][["beta"]][["coeffs"]]))
  }
})

test_that("every time varying cointegration algorithm reads the ML transition", {
  for (error in c("wishart", "gamma", "sv")) {

    set.seed(456789)
    plain <- vec_tvp_fitted(error, list(rho = 0.999))

    # A weight that leaves the transition at the identity in every direction is
    # no prior on the direction, and the chain is the one without it.
    set.seed(456789)
    negligible <- vec_tvp_fitted(error, list(rho = 0.999, p_tau_i = "ml", weight = 1e-12))
    expect_identical(negligible[["posterior"]][["beta"]][["coeffs"]],
                     plain[["posterior"]][["beta"]][["coeffs"]])

    # One that does not reaches the sampler through priors$beta$p_tau.
    set.seed(456789)
    informative <- suppressWarnings(
      vec_tvp_fitted(error, list(rho = 0.999, p_tau_i = "ml", weight = 0.1)))
    expect_false(is.null(informative[["priors"]][["beta"]][["p_tau"]]))
    expect_false(isTRUE(all.equal(informative[["posterior"]][["beta"]][["coeffs"]],
                                  plain[["posterior"]][["beta"]][["coeffs"]])))
  }
})
