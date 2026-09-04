test_that("VAR posterior draws have one row per stored iteration", {
  model <- fx_var_fitted()
  spec <- model[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_named(model[["posterior"]], c("a", "u_sigma_inv", "loglik"))
  expect_s3_class(model[["posterior"]][["a"]][["coeffs"]], "mcmc")
  expect_identical(dim(model[["posterior"]][["a"]][["coeffs"]]),
                   c(fx_iterations, n_coeffs))
  expect_equal(dim(model[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
               c(fx_iterations, spec[["k"]]^2))
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
})

test_that("the burn-in draws are discarded", {
  short <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = 12, burnin = 50)
  short <- add_priors(short, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 1, scale = 0.0001))
  short <- add_initial_values(short)
  set.seed(11)
  short <- add_posterior_coefficients(short)

  # Only the post-burn-in draws are kept, regardless of the burn-in length.
  expect_identical(nrow(short[["posterior"]][["a"]][["coeffs"]]), 12L)
})

test_that("posterior simulation is reproducible from the seed", {
  run <- function(seed) {
    set.seed(seed)
    add_posterior_coefficients(fx_var_initial())[["posterior"]][["a"]][["coeffs"]]
  }

  expect_equal(run(2024), run(2024))
  expect_false(isTRUE(all.equal(run(2024), run(2025))))
})

test_that("draws of the error precision are symmetric and positive definite", {
  draws <- fx_var_fitted()[["posterior"]][["u_sigma_inv"]][["coeffs"]]
  k <- fx_var_fitted()[["model"]][["k"]]

  for (i in seq_len(nrow(draws))) {
    sigma_inv <- matrix(draws[i, ], k)
    expect_equal(sigma_inv, t(sigma_inv))
    expect_true(all(eigen(sigma_inv, only.values = TRUE)[["values"]] > 0))
  }
})

test_that("the log-likelihood is stored per draw and observation", {
  model <- fx_var_fitted()
  loglik <- model[["posterior"]][["loglik"]]

  expect_s3_class(loglik, "mcmc")
  expect_identical(dim(loglik),
                   c(fx_iterations, nrow(model[["data"]][["train"]][["y"]])))
  expect_true(all(is.finite(loglik)))
})

test_that("the stored log-likelihood matches the normal density", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]
  y <- t(model[["data"]][["train"]][["y"]])
  x <- t(model[["data"]][["train"]][["x"]])

  draw <- 1L
  a <- matrix(model[["posterior"]][["a"]][["coeffs"]][draw, ], k)
  sigma_inv <- matrix(model[["posterior"]][["u_sigma_inv"]][["coeffs"]][draw, ], k)
  u <- y - a %*% x

  expected <- -0.5 * (k * log(2 * pi) - log(det(sigma_inv)) +
                        colSums(u * (sigma_inv %*% u)))
  expect_equal(as.numeric(model[["posterior"]][["loglik"]][draw, ]),
               as.numeric(expected))
})

test_that("VEC posterior draws cover alpha, beta and the error precision", {
  model <- fx_vec_fitted()
  spec <- model[["model"]]

  expect_true(all(c("a", "beta", "u_sigma_inv") %in% names(model[["posterior"]])))
  expect_identical(nrow(model[["posterior"]][["beta"]][["coeffs"]]),
                   fx_iterations)
  expect_identical(ncol(model[["posterior"]][["beta"]][["coeffs"]]),
                   as.integer(spec[["k_beta"]] * spec[["rank"]]))
  expect_true(all(is.finite(model[["posterior"]][["beta"]][["coeffs"]])))
})

test_that("posterior simulation runs for every model of a modellist", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                       sigma = list(df = 1, scale = 0.0001))
  models <- add_initial_values(models)
  set.seed(7)
  models <- add_posterior_coefficients(models)

  expect_s3_class(models, "modellist")
  expect_identical(
    vapply(models, function(x) ncol(x[["posterior"]][["a"]][["coeffs"]]),
           integer(1)),
    c(12L, 21L)
  )
})

test_that("stochastic volatility specifications produce draws", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            error = "sv", iterations = 10, burnin = 5)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1),
                      sigma = list(mu = 0, v_i = 0.01, shape = 3,
                                   rate = 0.0001, state_variance = 0.05,
                                   offset = 0.0001))
  model <- add_initial_values(model)
  set.seed(3)
  model <- add_posterior_coefficients(model)

  expect_identical(nrow(model[["posterior"]][["a"]][["coeffs"]]), 10L)
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
  # The volatility path is reported as one error precision matrix per period.
  tt <- nrow(model[["data"]][["train"]][["y"]])
  k <- model[["model"]][["k"]]
  expect_identical(dim(model[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
                   c(10L, as.integer(k * k * tt)))
})

test_that("a failing sampler is reported instead of aborting", {
  model <- fx_var_initial()
  # Corrupt the initial values so the sampler cannot run.
  model[["initial"]][["u_sigma_inv"]] <- matrix(NA_real_, 2, 2)

  # try() reports the sampler error on stderr, which is not part of the contract.
  utils::capture.output(
    failed <- add_posterior_coefficients(model),
    type = "message"
  )

  expect_true(isTRUE(failed[["error"]]))
  expect_null(failed[["posterior"]])
  # The object keeps its class so downstream code can filter on the flag.
  expect_s3_class(failed, "bvarmodel")
})

test_that("time varying parameters are drawn for every period", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            tvp = TRUE, iterations = 10, burnin = 5)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1,
                                  shape = 3, rate = 0.0001),
                      sigma = list(df = 3, scale = 0.0001))
  model <- add_initial_values(model)
  set.seed(3)
  model <- add_posterior_coefficients(model)

  spec <- model[["model"]]
  nobs <- nrow(model[["data"]][["train"]][["y"]])
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_identical(nrow(model[["posterior"]][["a"]][["coeffs"]]), 10L)
  expect_identical(ncol(model[["posterior"]][["a"]][["coeffs"]]),
                   as.integer(n_coeffs * nobs))
})
