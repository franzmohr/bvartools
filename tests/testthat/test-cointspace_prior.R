# cointspace_prior() builds the prior on the cointegration space for
# add_priors(), and for packages whose VEC models have a layout of their own.
# The priors it gives are tested with add_priors() in test-add_priors.R; what is
# tested here is what the exported function adds to that.

vec_coint_priors <- function(coint, model = fx_vec_model()) {
  add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10), coint = coint,
             sigma = list(df = "k", scale = 1))
}

test_that("add_priors() stores the prior cointspace_prior() builds", {
  model <- fx_vec_model()
  coints <- list(list(v_i = 0, p_tau_i = 1),
                 list(v_i = 0.1, p_tau_i = c(1, 2)),
                 list(v_i = 0.01, p_tau_i = "ml"),
                 list(v_i = "ml", p_tau_i = "ml", weight = 10))
  for (coint in coints) {
    expect_identical(vec_coint_priors(coint, model)[["priors"]][["beta"]],
                     cointspace_prior(model, coint))
  }
})

test_that("a model without cointegration gets no prior, but its coint is checked", {
  model <- fx_vec_model()
  model[["model"]][["rank"]] <- 0L
  expect_null(cointspace_prior(model, list(v_i = 0, p_tau_i = 1)))
  expect_error(cointspace_prior(model, list(v_i = 0)), "coint[$]p_tau_i")
  expect_error(cointspace_prior(model, list(v_i = 0, p_tau_i = 1, tau = 1)), "not recognised")
})

test_that("the prior is sized from the error correction term, not from model$m", {
  # A model of another layout, such as a sub-model of a global VEC model, whose
  # error correction term holds series that its model$m counts differently from
  # a model of create_bvecmodel().
  model <- fx_vec_model()
  k_w <- NCOL(model[["data"]][["train"]][["w"]])
  model[["model"]][["m"]] <- model[["model"]][["m"]] + 3L

  numeric_prior <- cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1))
  expect_identical(dim(numeric_prior[["p_tau_inv"]]), c(k_w, k_w))

  ml_prior <- cointspace_prior(model, list(v_i = 0.01, p_tau_i = "ml"))
  expect_identical(dim(ml_prior[["p_tau_inv"]]), c(k_w, k_w))

  expect_error(cointspace_prior(model, list(v_i = 0.1, p_tau_i = diag(k_w + 3))),
               paste0(k_w, " x ", k_w))
})

test_that("a cointegration rank without an error correction term is refused", {
  model <- fx_vec_model()
  model[["data"]][["train"]][["w"]] <- NULL
  expect_error(cointspace_prior(model, list(v_i = 0, p_tau_i = 1)), "error correction term")
})

# coint$g_i: the fixed G^-1 that VecNormalStochvol scales its loadings' prior by.

vec_sv_model <- function(tvp = FALSE, error = "sv") {
  create_bvecmodel(vec_data(), p = 1, r = 1, const = "unrestricted",
                   tvp = tvp, error = error,
                   iterations = fx_iterations, burnin = fx_burnin)
}

vec_sv_fitted <- function(coint, seed = 20260921) {
  object <- add_priors(vec_sv_model(),
                       coef = list(v_i = 1, v_i_det = 1 / 10),
                       coint = coint,
                       sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                    state_variance = 0.05, offset = 1e-8))
  object <- add_initial_values(object)
  add_posterior_coefficients(add_seed(object, seed))
}

test_that("coint$g_i is built as a k x k precision", {
  model <- vec_sv_model()
  k <- NCOL(model[["data"]][["train"]][["y"]])

  expect_null(cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1))[["g_inv"]])
  expect_identical(cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1, g_i = 2))[["g_inv"]],
                   diag(2, k))
  expect_identical(cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1,
                                                g_i = seq_len(k)))[["g_inv"]],
                   diag(as.numeric(seq_len(k)), k))
  full <- matrix(0.5, k, k) + diag(1, k)
  expect_identical(cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1, g_i = full))[["g_inv"]],
                   full)

  # "ml" is the inverse of Johansen's error covariance, the one v_i = "ml" uses.
  ml <- bvartools:::.coint_ml(model)
  expect_equal(cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1, g_i = "ml"))[["g_inv"]],
               solve(ml[["omega"]]))

  # sv+covar has a G of its own to give as well.
  expect_identical(cointspace_prior(vec_sv_model(error = "sv+covar"),
                                    list(v_i = 0.1, p_tau_i = 1, g_i = 1))[["g_inv"]],
                   diag(1, k))
})

test_that("coint$g_i is refused wherever G is not fixed by the file", {
  coint <- list(v_i = 0.1, p_tau_i = 1, g_i = 1)
  expect_error(cointspace_prior(fx_vec_model(), coint), "stochastic volatility")
  expect_error(cointspace_prior(vec_sv_model(error = "gamma"), coint), "stochastic volatility")
  expect_error(cointspace_prior(vec_sv_model(tvp = TRUE), c(coint, rho = 0.999)),
               "stochastic volatility")
})

test_that("a coint$g_i that is not a precision is refused", {
  model <- vec_sv_model()
  k <- NCOL(model[["data"]][["train"]][["y"]])
  prior <- function(g_i) cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1, g_i = g_i))

  expect_error(prior("identity"), "numeric or \"ml\"")
  expect_error(prior(rep(1, k + 1)), "one per endogenous variable")
  expect_error(prior(diag(k + 1)), paste0(k, " x ", k))
  asymmetric <- diag(k)
  asymmetric[1, 2] <- 0.5
  expect_error(prior(asymmetric), "symmetric")
  expect_error(prior(0), "positive definite")
  expect_error(prior(-1), "positive definite")
  expect_error(prior(matrix(1, k, k)), "positive definite")
})

test_that("coint$g_i reaches the sampler", {
  without <- vec_sv_fitted(list(v_i = 0.1, p_tau_i = 1))
  with_g <- vec_sv_fitted(list(v_i = 0.1, p_tau_i = 1, g_i = 1e4))
  expect_equal(with_g[["model"]][["algorithm"]], "VecNormalStochvol")
  expect_equal(with_g[["priors"]][["beta"]][["g_inv"]], diag(1e4, 2))

  # Same seed, same starting values: only G differs, and a G^-1 this large
  # shrinks the loadings hard towards zero.
  expect_false(isTRUE(all.equal(as.matrix(with_g[["posterior"]][["a"]][["coeffs"]]),
                                as.matrix(without[["posterior"]][["a"]][["coeffs"]]))))
  expect_lt(mean(abs(as.matrix(with_g[["posterior"]][["a"]][["coeffs"]])[, 1:2])),
            mean(abs(as.matrix(without[["posterior"]][["a"]][["coeffs"]])[, 1:2])))

  # And the binding hands it to the core rather than dropping it: a g_inv put
  # by hand on a model that must not have one is refused by validate().
  object <- add_initial_values(fx_vec_priors())
  object[["priors"]][["beta"]][["g_inv"]] <- diag(2)
  expect_error(add_posterior_coefficients(object), "VecNormalStochvol alone")
})

test_that("coint$g_i survives the HDF5 round trip under the name BayesTS reads", {
  skip_if_not_installed("hdf5r")
  object <- add_priors(vec_sv_model(),
                       coef = list(v_i = 1, v_i_det = 1 / 10),
                       coint = list(v_i = 0.1, p_tau_i = 1, g_i = c(2, 3)),
                       sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                    state_variance = 0.05, offset = 1e-8))
  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)

  file <- hdf5r::H5File$new(path, mode = "r")
  expect_true(file$exists("priors/beta/g_inv"))
  file$close_all()

  restored <- read_model_from_hdf5(path)
  expect_equal(restored[["priors"]][["beta"]][["g_inv"]], diag(c(2, 3)),
               ignore_attr = TRUE)
})
