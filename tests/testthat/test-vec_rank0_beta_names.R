# A VEC model of rank zero needs no prior on a cointegration space, and the draws
# of the cointegration vectors are named after what they weight.

rank_models <- function(r, ...) {
  create_bvecmodel(vec_data(), p = 2, r = r, const = "unrestricted",
                   iterations = fx_iterations, burnin = fx_burnin, ...)
}
flat <- list(v_i = 0, v_i_det = 0)
wishart <- list(df = "k", scale = 1)
uniform <- list(v_i = 0, p_tau_i = 1)

test_that("a VEC model of rank zero takes no cointegration prior", {
  without <- add_priors(rank_models(0), coef = flat, sigma = wishart)
  with <- add_priors(rank_models(0), coef = flat, coint = uniform, sigma = wishart)
  expect_null(without[["priors"]][["beta"]])
  expect_identical(without[["priors"]], with[["priors"]])
  fitted <- add_posterior_coefficients(add_initial_values(without))
  expect_null(fitted[["posterior"]][["beta"]])

  # A prior that is given is still checked.
  expect_error(add_priors(rank_models(0), coef = flat, coint = list(nonsense = 1), sigma = wishart),
               "not recognised")
})

test_that("a positive rank still needs one, and says so", {
  expect_error(add_priors(rank_models(1), coef = flat, sigma = wishart),
               "must be specified for a VEC model with a cointegration rank of 1")
})

test_that("a list of models over ranks shares one set of priors", {
  models <- add_priors(rank_models(0:2), coef = flat, coint = uniform, sigma = wishart)
  expect_null(models[[1]][["priors"]][["beta"]])
  expect_false(is.null(models[[2]][["priors"]][["beta"]]))
  expect_error(add_priors(rank_models(0:2), coef = flat, sigma = wishart), "rank of 1")
})

test_that("the draws of the cointegration vectors are named after the series they weight", {
  fitted <- add_posterior_coefficients(add_seed(add_initial_values(
    add_priors(rank_models(2), coef = flat, coint = uniform, sigma = wishart)), 31))
  w <- colnames(fitted[["data"]][["train"]][["w"]])
  expect_identical(colnames(fitted[["posterior"]][["beta"]][["coeffs"]]),
                   c(paste0("ect1.", w), paste0("ect2.", w)))

  # They travel through a file and through several chains.
  path <- temp_h5_file()
  write_to_hdf5(fitted, path)
  back <- read_model_from_hdf5(path)
  expect_identical(colnames(back[["posterior"]][["beta"]][["coeffs"]]),
                   colnames(fitted[["posterior"]][["beta"]][["coeffs"]]))
  expect_equal(back[["posterior"]][["beta"]][["coeffs"]], fitted[["posterior"]][["beta"]][["coeffs"]])

  two <- add_posterior_coefficients(add_seed(add_initial_values(
    add_priors(rank_models(1), coef = flat, coint = uniform, sigma = wishart)), 31), chains = 2)
  expect_identical(colnames(two[["posterior"]][["beta"]][["coeffs"]]), paste0("ect1.", w))
})

test_that("time varying cointegration vectors are named by period as well", {
  fitted <- fx_vec_tvp_fitted("gamma")
  beta <- fitted[["posterior"]][["beta"]][["coeffs"]]
  w <- colnames(fitted[["data"]][["train"]][["w"]])
  tt <- nrow(fitted[["data"]][["train"]][["y"]])
  expect_identical(ncol(beta), tt * length(w))
  expect_identical(colnames(beta)[seq_along(w)], paste0("ect1.", w, ".t1"))
  expect_identical(utils::tail(colnames(beta), 1), paste0("ect1.", utils::tail(w, 1), ".t", tt))
})
