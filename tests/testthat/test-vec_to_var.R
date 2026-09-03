test_that("vec_to_var produces a VAR model of the same size", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)

  expect_s3_class(var, "bvarmodel")
  expect_identical(var[["model"]][["type"]], "VAR")
  expect_identical(var[["model"]][["k"]], vec[["model"]][["k"]])
  expect_identical(var[["model"]][["endogen"]], vec[["model"]][["endogen"]])
  # p is the lag order of the model in levels in both parameterisations.
  expect_equal(var[["model"]][["p"]], vec[["model"]][["p"]])
})

test_that("the converted model carries the posterior draws over", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)
  spec <- var[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_identical(nrow(var[["posterior"]][["a"]][["coeffs"]]), fx_iterations)
  expect_identical(ncol(var[["posterior"]][["a"]][["coeffs"]]),
                   as.integer(n_coeffs))
  # The error precision is unchanged by the reparameterisation.
  expect_equal(var[["posterior"]][["u_sigma_inv"]][["coeffs"]],
               vec[["posterior"]][["u_sigma_inv"]][["coeffs"]])
})

test_that("the level coefficients follow from Pi and the short-run terms", {
  vec <- fx_vec_fitted()
  var <- vec_to_var(vec)
  k <- vec[["model"]][["k"]]
  r <- vec[["model"]][["rank"]]

  draw <- 1L
  alpha <- matrix(vec[["posterior"]][["a"]][["coeffs"]][draw, seq_len(k * r)], k)
  beta <- matrix(vec[["posterior"]][["beta"]][["coeffs"]][draw, ],
                 vec[["model"]][["k_beta"]])

  # With p = 1 there are no lagged differences, so A_1 = Pi + I.
  expected <- alpha %*% t(beta) + diag(1, k)
  actual <- matrix(var[["posterior"]][["a"]][["coeffs"]][draw, seq_len(k^2)], k)
  expect_equal(actual, expected)
})

test_that("the converted model is usable by the application functions", {
  var <- vec_to_var(fx_vec_fitted())

  expect_s3_class(irf(var, impulse = "R", response = "Dp", n_ahead = 3),
                  "bvarirf")
  expect_s3_class(fevd(var, response = "Dp", n_ahead = 3), "bvarfevd")
  expect_no_error(summary(var))
})

test_that("a modellist of VEC models is converted element by element", {
  models <- create_bvecmodel(vec_data(), p = 1, r = 1:2,
                             const = "unrestricted",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 1, v_i_det = 1 / 10),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = "k", scale = 1))
  models <- add_initial_values(models)
  set.seed(17)
  models <- add_posterior_coefficients(models)
  converted <- vec_to_var(models)

  expect_s3_class(converted, "modellist")
  expect_length(converted, 2)
  expect_true(all(vapply(converted, inherits, logical(1), "bvarmodel")))
})
