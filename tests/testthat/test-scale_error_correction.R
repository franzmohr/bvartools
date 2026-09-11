# The Koop et al. (2010) prior lives on the cointegration space, so the scale of
# the series in the error correction term matters. scale_error_correction()
# divides them by the standard deviation of the corresponding differences and
# rescale_error_correction() undoes that once the draws are in.

# A fitted model whose error correction term was scaled before estimation.
scaled_vec_fitted <- function(const = "unrestricted") {
  model <- create_bvecmodel(vec_data(), p = 1, r = 1, const = const,
                            iterations = fx_iterations, burnin = fx_burnin)
  model <- scale_error_correction(model)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1))
  model <- add_initial_values(model)
  set.seed(456789)
  add_posterior_coefficients(model)
}

test_that("scaling divides the error correction term by the scaling factors", {
  model <- scale_error_correction(fx_vec_model())
  raw <- fx_vec_model()[["data"]][["train"]][["w"]]
  factors <- attr(model[["data"]][["train"]][["w"]], "scale")

  expect_identical(names(factors), dimnames(raw)[[2]])
  expect_equal(unname(factors), unname(apply(diff(raw), 2, stats::sd)))
  expect_equal(as.numeric(model[["data"]][["train"]][["w"]]),
               as.numeric(raw / t(matrix(factors, length(factors), nrow(raw)))))
})

test_that("rescaling restores the error correction term", {
  fitted <- scaled_vec_fitted()
  rescaled <- rescale_error_correction(fitted)

  expect_equal(as.numeric(rescaled[["data"]][["train"]][["w"]]),
               as.numeric(fx_vec_model()[["data"]][["train"]][["w"]]))
  # Nothing is left to undo, so a second call is an error rather than a second
  # transformation.
  expect_null(attr(rescaled[["data"]][["train"]][["w"]], "scale"))
  expect_error(rescale_error_correction(rescaled), "attribute 'scale'")
})

test_that("rescaling transforms beta and leaves alpha and Pi alone", {
  fitted <- scaled_vec_fitted()
  rescaled <- rescale_error_correction(fitted)
  k <- fitted[["model"]][["k"]]
  r <- fitted[["model"]][["rank"]]
  k_beta <- fitted[["model"]][["k_beta"]]
  factors <- attr(fitted[["data"]][["train"]][["w"]], "scale")

  alpha <- fitted[["posterior"]][["a"]][["coeffs"]][, seq_len(k * r), drop = FALSE]
  # The loadings do not depend on the scale of the error correction term.
  expect_equal(rescaled[["posterior"]][["a"]][["coeffs"]], fitted[["posterior"]][["a"]][["coeffs"]])

  draw <- 1L
  beta_scaled <- matrix(fitted[["posterior"]][["beta"]][["coeffs"]][draw, ], k_beta)
  beta_level <- matrix(rescaled[["posterior"]][["beta"]][["coeffs"]][draw, ], k_beta)
  expect_equal(beta_level, diag(1 / factors, k_beta) %*% beta_scaled)

  # Pi %*% w is the same number in both parameterisations.
  alpha_draw <- matrix(alpha[draw, ], k)
  w_scaled <- fitted[["data"]][["train"]][["w"]]
  w_level <- rescaled[["data"]][["train"]][["w"]]
  expect_equal(w_level %*% beta_level %*% t(alpha_draw),
               w_scaled %*% beta_scaled %*% t(alpha_draw))
})

test_that("rescaling works with deterministic terms in the cointegration space", {
  fitted <- scaled_vec_fitted(const = "restricted")
  rescaled <- rescale_error_correction(fitted)

  expect_identical(dimnames(rescaled[["data"]][["train"]][["w"]])[[2]],
                   c("l.R", "l.Dp", "const"))
  # The restricted constant is a column of ones again.
  expect_equal(as.numeric(rescaled[["data"]][["train"]][["w"]][, "const"]),
               rep(1, nrow(rescaled[["data"]][["train"]][["w"]])))
})

test_that("vec_to_var refuses a model whose error correction term is scaled", {
  fitted <- scaled_vec_fitted()

  expect_error(vec_to_var(fitted), "rescale_error_correction")
  expect_no_error(vec_to_var(rescale_error_correction(fitted)))
})

test_that("vec_to_var recovers the levels of the endogenous variables", {
  var <- vec_to_var(rescale_error_correction(scaled_vec_fitted()))
  levels <- vec_data()

  expect_equal(as.numeric(var[["data"]][["train"]][["y"]]),
               as.numeric(stats::window(levels, start = stats::start(var[["data"]][["train"]][["y"]]))))
})

test_that("bvec reconstructs the levels from a scaled error correction term", {
  fitted <- scaled_vec_fitted()
  k <- fitted[["model"]][["k"]]
  r <- fitted[["model"]][["rank"]]
  draws <- t(fitted[["posterior"]][["a"]][["coeffs"]])
  n_a <- k * r

  object <- bvec(y = fitted[["data"]][["train"]][["y"]],
                 w = fitted[["data"]][["train"]][["w"]],
                 x_d = fitted[["data"]][["train"]][["x"]],
                 r = r,
                 alpha = draws[seq_len(n_a), , drop = FALSE],
                 beta = t(fitted[["posterior"]][["beta"]][["coeffs"]]),
                 C = draws[-seq_len(n_a), , drop = FALSE],
                 Sigma = t(.invert_sigma_draws(t(fitted[["posterior"]][["u_sigma_inv"]][["coeffs"]]), k)))

  # The levels are the differences plus their first lag, which only holds if the
  # scaling of the error correction term is undone first.
  expect_equal(as.numeric(object[["data"]][["original"]][["endogen"]]),
               as.numeric(stats::window(vec_data(),
                                        start = stats::start(object[["data"]][["train"]][["y"]]))))
})

test_that("a second scaling is refused", {
  # Scaling twice does not change the numbers -- dividing by the standard
  # deviation of the differences makes that standard deviation one, so the
  # second pass divides by one -- but it recomputes the factors from the series
  # it is given and overwrites the stored ones, which makes the first scaling
  # irreversible.
  model <- scale_error_correction(fx_vec_model())

  expect_error(scale_error_correction(model), "are already scaled")
})

test_that("rescaling does not require posterior draws", {
  # A model is scaled before it is estimated when it is exported for an
  # external sampler, and a sub-model whose run produced nothing is left in the
  # same state. Neither could be put back when rescaling insisted on draws.
  model <- fx_vec_model()
  scaled <- scale_error_correction(model)

  expect_null(scaled[["posterior"]])

  restored <- rescale_error_correction(scaled)

  expect_null(attr(restored[["data"]][["train"]][["w"]], "scale"))
  expect_equal(unclass(restored[["data"]][["train"]][["w"]]),
               unclass(model[["data"]][["train"]][["w"]]),
               ignore_attr = TRUE)
})

test_that("the scaling factors survive a write and read round trip", {
  skip_if_not_installed("hdf5r")

  # Without them an exported model that was scaled could never be put back on
  # the scale of the data, so the draws an external sampler writes into it
  # would not be interpretable.
  model <- scale_error_correction(fx_vec_model())
  factors <- attr(model[["data"]][["train"]][["w"]], "scale")

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(attr(restored[["data"]][["train"]][["w"]], "scale"), factors)
  expect_equal(unclass(restored[["data"]][["train"]][["w"]]),
               unclass(model[["data"]][["train"]][["w"]]),
               ignore_attr = TRUE)

  # And what came back can be put back on the scale of the data.
  expect_equal(unclass(rescale_error_correction(restored)[["data"]][["train"]][["w"]]),
               unclass(rescale_error_correction(model)[["data"]][["train"]][["w"]]),
               ignore_attr = TRUE)
})

test_that("an unscaled model gets no scaling factors written", {
  skip_if_not_installed("hdf5r")

  path <- temp_h5_file()
  write_to_hdf5(fx_vec_fitted(), filename = path)
  restored <- read_model_from_hdf5(path)

  expect_null(attr(restored[["data"]][["train"]][["w"]], "scale"))
})
