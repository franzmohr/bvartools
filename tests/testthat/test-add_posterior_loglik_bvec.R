# The log-likelihood of an error correction model.
#
# Every algorithm a 'bvecmodel' can be estimated with has a log-likelihood
# implemented in C++. The method used to admit only one of them, which left the
# selection criteria unavailable for the others and made a specification search
# over error and coefficient processes impossible.

vec_fitted <- function(tvp = FALSE, error = "wishart", iterations = 30,
                       burnin = 15) {

  data("us_macrodata", envir = environment())

  object <- create_bvecmodel(data = us_macrodata, p = 2,
                             const = "unrestricted", r = 1,
                             tvp = tvp, error = error,
                             iterations = iterations, burnin = burnin)

  coef_prior <- list(v_i = 1, v_i_det = 1 / 10)
  if (tvp) {
    coef_prior <- c(coef_prior, list(shape = 3, rate = 1e-8, rate_det = 1e-8))
  }
  sigma_prior <- if (error %in% c("sv", "sv+covar")) {
    list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
         state_variance = 0.05, offset = 1e-8)
  } else if (error %in% c("gamma", "gamma+covar")) {
    list(shape = 3, rate = 0.01)
  } else {
    list(df = 3, scale = 1)
  }

  object <- add_priors(object, coef = coef_prior, sigma = sigma_prior,
                       coint = if (tvp) list(rho = 0.999)
                               else list(v_i = 0, p_tau_i = 1))
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

test_that("every error correction algorithm has a log-likelihood", {
  specifications <- list(
    list(tvp = FALSE, error = "wishart",   algorithm = "VecNormalWishart"),
    list(tvp = FALSE, error = "sv",        algorithm = "VecNormalStochvol"),
    list(tvp = FALSE, error = "gamma",     algorithm = "VecNormalGamma"),
    list(tvp = TRUE,  error = "wishart",   algorithm = "VecTvpWishart"),
    list(tvp = TRUE,  error = "sv",        algorithm = "VecTvpStochvol"),
    list(tvp = TRUE,  error = "gamma",     algorithm = "VecTvpGamma")
  )

  for (specification in specifications) {

    object <- vec_fitted(tvp = specification[["tvp"]],
                         error = specification[["error"]])

    expect_equal(object[["model"]][["algorithm"]],
                 specification[["algorithm"]])

    object <- add_posterior_loglik(object)
    loglik <- object[["posterior"]][["loglik"]]

    # One value per draw and period, and every one of them usable: a criterion
    # built on an infinite or missing log-likelihood ranks nothing.
    expect_equal(nrow(loglik),
                 nrow(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
                 info = specification[["algorithm"]])
    expect_equal(ncol(loglik),
                 nrow(object[["data"]][["train"]][["y"]]),
                 info = specification[["algorithm"]])
    expect_true(all(is.finite(as.matrix(loglik))),
                info = specification[["algorithm"]])

    # Which is what the selection criteria are for.
    criteria <- selection_criteria(object)
    expect_true(is.finite(criteria[["LL"]][["median"]]),
                info = specification[["algorithm"]])
    expect_true(is.finite(criteria[["AIC"]][["median"]]),
                info = specification[["algorithm"]])
  }
})

test_that("an algorithm without a log-likelihood is refused by name", {
  object <- vec_fitted()
  object[["model"]][["algorithm"]] <- "SomethingElse"

  expect_error(add_posterior_loglik(object), "Algorithm not implemented yet")
})

test_that("the log-likelihood needs the draws of the precision", {
  object <- vec_fitted()
  object[["posterior"]][["u_sigma_inv"]][["coeffs"]] <- NULL

  expect_error(add_posterior_loglik(object),
               "does not contain posterior draws")
})
