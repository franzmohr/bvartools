# expected_model_size() and the warnings built on it.
#
# The size is predicted from the specification, so the invariant is that the
# prediction is what the samplers then store: for every algorithm and every
# option that adds a block, each element expected_model_size() lists is in the
# fitted posterior with exactly the predicted dimensions, and nothing in the
# posterior is missing from the list. A sampler that starts storing something
# new fails here. The warnings of add_posterior_coefficients() and
# write_to_hdf5() are checked for being raised once per call, not per model.

size_iterations <- 7L

size_sigma_prior <- function(error) {
  switch(sub("\\+covar", "", error),
         sv = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                   state_variance = 0.05, offset = 1e-8),
         gamma = list(shape = 3, rate = 0.01),
         ald = list(shape = 3, rate = 0.01),
         list(df = 3, scale = 1))
}

size_var <- function(error = "wishart", tvp = FALSE, varsel = "none", structural = FALSE,
                     coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                     sigma = size_sigma_prior(error), ...) {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", error = error,
                            tvp = tvp, varsel = varsel, structural = structural,
                            iterations = size_iterations, burnin = 2, ...)
  varsel_prior <- NULL
  if (varsel != "none") {
    varsel_prior <- list(inprior = 0.5, covar = grepl("covar", error))
    if (varsel == "ssvs") {
      varsel_prior[["semiautomatic"]] <- c(0.1, 10)
    }
  }
  suppressWarnings(add_initial_values(add_priors(model, coef = coef, sigma = sigma,
                                                 varsel = varsel_prior)))
}

size_vec <- function(error = "wishart", tvp = FALSE, r = 1, coint = NULL,
                     algorithm = NULL, varsel = "none") {
  model <- create_bvecmodel(stats::window(at_domestic()[, c("lr", "Dp", "r")], end = c(2005, 4)) * 100,
                            p = 2, r = r, const = "unrestricted", error = error, tvp = tvp,
                            algorithm = algorithm, varsel = varsel,
                            iterations = size_iterations, burnin = 2)
  if (is.null(coint)) {
    coint <- if (tvp) list(rho = 0.999) else list(v_i = 0, p_tau_i = 1)
  }
  suppressWarnings(add_initial_values(add_priors(
    model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-8, rate_det = 1e-8),
    sigma = size_sigma_prior(error), coint = coint,
    varsel = if (varsel != "none") list(inprior = 0.5, covar = FALSE))))
}

# The dimensions of every matrix in the posterior, as "element draws columns".
posterior_dims <- function(object) {
  post <- object[["posterior"]]
  result <- character(0)
  for (block in names(post)) {
    if (is.list(post[[block]])) {
      for (element in names(post[[block]])) {
        x <- post[[block]][[element]]
        result <- c(result, paste(paste0("posterior$", block, "$", element), NROW(x), NCOL(x)))
      }
    } else {
      result <- c(result, paste(paste0("posterior$", block), NROW(post[[block]]), NCOL(post[[block]])))
    }
  }
  sort(result)
}

predicted_dims <- function(size, step = "add_posterior_coefficients") {
  size <- size[size[["step"]] == step, ]
  sort(paste(size[["element"]], size[["draws"]], size[["columns"]]))
}

expect_size_matches <- function(model, label) {
  predicted <- expected_model_size(model)
  set.seed(1)
  fitted <- suppressWarnings(add_posterior_coefficients(model))
  expect_identical(predicted_dims(predicted), posterior_dims(fitted), info = label)
}

test_that("the predicted blocks are what every VAR sampler stores", {
  for (tvp in c(FALSE, TRUE)) {
    for (error in c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      expect_size_matches(size_var(error, tvp), paste("VAR", error, tvp))
    }
    expect_size_matches(size_var("ald", tvp, quantile = 0.5), paste("VAR ald", tvp))
    for (error in c("wishart", "gamma+covar", "sv+covar")) {
      expect_size_matches(size_var(error, tvp, varsel = "bvs"), paste("VAR bvs", error, tvp))
    }
    expect_size_matches(size_var("gamma", tvp, structural = TRUE), paste("SVAR", tvp))
  }
  expect_size_matches(size_var(varsel = "ssvs"), "VAR ssvs")
  expect_size_matches(
    size_var("sv+covar", TRUE, coef = list(v_i = 1, v_i_det = 0.1, omega_v = 0.001),
             sigma = list(mu = 0, v_i = 0.01, omega_v = 0.1, state_variance = 0.05, offset = 1e-8)),
    "VAR non-centred")
})

test_that("the predicted blocks are what every VEC sampler stores", {
  for (tvp in c(FALSE, TRUE)) {
    for (error in c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      expect_size_matches(size_vec(error, tvp), paste("VEC", error, tvp))
    }
    expect_size_matches(size_vec(tvp = tvp, varsel = "bvs"), paste("VEC bvs", tvp))
  }
  expect_size_matches(size_vec(r = 2), "VEC rank 2")
  expect_size_matches(size_vec(r = 0), "VEC rank 0")
  expect_size_matches(size_vec("sv", TRUE, r = 0), "VEC rank 0 sv tvp")
  expect_size_matches(size_vec(algorithm = "KLGS2010"), "VEC KLGS2010")
  expect_size_matches(size_vec("sv", TRUE, coint = list(rho = 0.99, rho_min = 0.9, rho_max = 0.999)),
                      "VEC drawn rho")
})

test_that("the discounted models are sized by period rather than by draw", {
  var <- create_bvarmodel(var_data(), p = 1, deterministic = "const", algorithm = "discount",
                          iterations = 10, burnin = 0)
  var <- add_initial_values(add_priors(var, coef = list(v_i = 1, v_i_det = 0.1),
                                       sigma = list(df = "k", scale = 1)))
  expect_size_matches(var, "VarTvpDiscount")

  vec <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                          algorithm = "discount", iterations = 10, burnin = 0)
  vec <- add_initial_values(add_priors(vec, coef = list(v_i = 1, v_i_det = 0.1),
                                       sigma = list(df = "k", scale = 1)))
  expect_size_matches(vec, "VecTvpDiscount")
})

test_that("the log-likelihood and the forecasts are predicted as well", {
  model <- add_forecast_input(size_var("sv", TRUE), n_ahead = 4)
  size <- expected_model_size(model)

  set.seed(2)
  fitted <- add_posterior_coefficients(model)
  fitted <- add_posterior_loglik(fitted)
  fitted <- add_posterior_forecasts(fitted)

  expect_identical(predicted_dims(size, "add_posterior_loglik"),
                   paste("posterior$loglik", paste(dim(fitted[["posterior"]][["loglik"]]), collapse = " ")))
  expect_identical(predicted_dims(size, "add_posterior_forecasts"),
                   paste("posterior$forecast$forecasts",
                         paste(dim(fitted[["posterior"]][["forecast"]][["forecasts"]]), collapse = " ")))

  # Eight bytes per number, and in total what the fitted object then takes, up
  # to the attributes each matrix carries, which do not grow with the draws.
  expect_equal(size[["bytes"]][-1], 8 * size[["draws"]][-1] * size[["columns"]][-1])
  actual <- as.numeric(utils::object.size(fitted))
  expect_lt(actual - sum(size[["bytes"]]), 1000 * nrow(size))
  expect_gte(actual, sum(size[["bytes"]]))
})

test_that("draws are counted once per chain", {
  model <- size_var()
  one <- expected_model_size(model)
  three <- expected_model_size(model, chains = 3)
  expect_equal(three[["draws"]][-1], 3 * one[["draws"]][-1])

  fitted <- add_posterior_coefficients(add_seed(model, 7), chains = 3)
  expect_equal(NROW(fitted[["posterior"]][["a"]][["coeffs"]]),
               three[["draws"]][three[["element"]] == "posterior$a$coeffs"])
})

test_that("a list of models is sized model by model", {
  models <- fx_var_modellist()
  size <- expected_model_size(models)
  expect_s3_class(size, "modelsize")
  expect_setequal(unique(size[["model"]]), c("1", "2"))
  expect_equal(sum(size[["bytes"]]),
               sum(expected_model_size(models[[1]])[["bytes"]]) + sum(expected_model_size(models[[2]])[["bytes"]]))
  expect_output(print(size), "2 models")
  expect_output(print(expected_model_size(models[[1]])), "posterior\\$a\\$coeffs")

  expect_error(expected_model_size(list(1)), "no method")
})

test_that("drawing a model that is too large warns once, before it starts", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = size_iterations, burnin = 2)
  models <- add_initial_values(add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                                          sigma = list(df = 1, scale = 0.0001)))

  old <- options(bvartools.size_warning = 1000)
  on.exit(options(old), add = TRUE)

  warnings <- character(0)
  withCallingHandlers(add_posterior_coefficients(models),
                      warning = function(w) {
                        warnings <<- c(warnings, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_length(warnings, 1)
  expect_match(warnings, "list of models")
  expect_match(warnings, "bvartools.size_warning = Inf")

  options(bvartools.size_warning = Inf)
  expect_no_warning(add_posterior_coefficients(models[[1]]))
})

test_that("writing an object that is too large warns once", {
  models <- fx_var_modellist()
  folder <- tempfile("size")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  old <- options(bvartools.size_warning = 1000)
  on.exit(options(old), add = TRUE)

  warnings <- character(0)
  withCallingHandlers(write_to_hdf5(models, folder),
                      warning = function(w) {
                        warnings <<- c(warnings, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_length(warnings, 1)
  expect_match(warnings, "disk space")
})
