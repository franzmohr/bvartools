# Centring the error correction term. scale_error_correction(centre = TRUE)
# subtracts the means of the stochastic series before the simulation, the
# unrestricted constant takes them up, and rescale_error_correction() writes the
# draws back in terms of the series as they are.

centred_vec_model <- function(tvp = FALSE, const = "unrestricted") {
  create_bvecmodel(vec_data(), p = 2, r = 1, const = const,
                   trend = if (tvp) NULL else "restricted", tvp = tvp,
                   error = if (tvp) "gamma" else "wishart",
                   iterations = fx_iterations, burnin = fx_burnin)
}

centred_vec_priors <- function(model, tvp = FALSE) {
  if (tvp) {
    add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
               coint = list(rho = 0.99), sigma = tvp_sigma_prior("gamma"))
  } else {
    add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = "k", scale = 1))
  }
}

centred_vec_fitted <- function(tvp = FALSE, scale = FALSE) {
  cached_fixture(paste0("vec_centred_", tvp, "_", scale), {
    model <- scale_error_correction(centred_vec_model(tvp), scale = scale, centre = TRUE)
    set.seed(20260914)
    model <- add_initial_values(centred_vec_priors(model, tvp))
    add_posterior_loglik(add_posterior_coefficients(model))
  })
}

# The fitted values of one draw, period by period: alpha beta' w + C x with the
# coefficients of each period, or of the only one.
fitted_of <- function(object, draw = 1) {
  k <- object[["model"]][["k"]]
  r <- object[["model"]][["rank"]]
  w <- object[["data"]][["train"]][["w"]]
  x <- object[["data"]][["train"]][["x"]]
  k_ect <- ncol(w)
  n_coef <- k * r + k * ncol(x)
  a <- object[["posterior"]][["a"]][["coeffs"]][draw, ]
  beta <- object[["posterior"]][["beta"]][["coeffs"]][draw, ]
  a_periods <- length(a) / n_coef
  beta_periods <- length(beta) / (k_ect * r)

  t(vapply(seq_len(nrow(w)), function(t) {
    a_t <- matrix(a[(min(t, a_periods) - 1) * n_coef + seq_len(n_coef)], k)
    beta_t <- matrix(beta[(min(t, beta_periods) - 1) * k_ect * r + seq_len(k_ect * r)], k_ect)
    as.numeric(a_t[, seq_len(r), drop = FALSE] %*% t(beta_t) %*% as.numeric(w[t, ]) +
                 a_t[, -seq_len(r), drop = FALSE] %*% as.numeric(x[t, ]))
  }, numeric(k)))
}

test_that("centring subtracts the means of the stochastic series and nothing else", {
  model <- centred_vec_model()
  raw <- model[["data"]][["train"]][["w"]]
  centred <- scale_error_correction(model, scale = FALSE, centre = TRUE)
  w <- centred[["data"]][["train"]][["w"]]
  means <- attr(w, "centre")

  expect_identical(names(means), dimnames(raw)[[2]])
  expect_null(attr(w, "scale"))
  # The restricted trend is a deterministic term and keeps its values.
  expect_equal(unname(means[["trend"]]), 0)
  expect_equal(as.numeric(w[, "trend"]), as.numeric(raw[, "trend"]))
  stochastic <- dimnames(raw)[[2]] != "trend"
  expect_equal(unname(means[stochastic]), unname(colMeans(raw[, stochastic])))
  expect_equal(unname(colMeans(w[, stochastic])), rep(0, sum(stochastic)))
  expect_identical(stats::tsp(w), stats::tsp(raw))
})

test_that("centring and scaling centre first and scale second", {
  model <- centred_vec_model()
  raw <- as.matrix(model[["data"]][["train"]][["w"]])
  both <- scale_error_correction(model, centre = TRUE)
  w <- both[["data"]][["train"]][["w"]]
  factors <- attr(w, "scale")
  means <- attr(w, "centre")

  expect_equal(unname(factors), unname(attr(scale_error_correction(model)[["data"]][["train"]][["w"]], "scale")))
  expect_equal(as.numeric(w), as.numeric(t((t(raw) - means) / factors)))
})

test_that("centring needs an unrestricted constant and a flag to act on", {
  expect_error(scale_error_correction(centred_vec_model(const = "restricted"), centre = TRUE),
               "unrestricted constant")
  expect_error(scale_error_correction(centred_vec_model(), scale = FALSE), "At least one")
  expect_error(scale_error_correction(centred_vec_model(), centre = NA), "TRUE or FALSE")
})

test_that("a second transformation is refused", {
  centred <- scale_error_correction(centred_vec_model(), scale = FALSE, centre = TRUE)

  expect_error(scale_error_correction(centred, scale = FALSE, centre = TRUE), "already scaled or centred")
  expect_error(scale_error_correction(centred), "already scaled or centred")
})

test_that("rescaling a centred model keeps the fitted values of every draw", {
  for (tvp in c(FALSE, TRUE)) {
    fitted <- centred_vec_fitted(tvp = tvp)
    back <- rescale_error_correction(fitted)
    draws <- nrow(fitted[["posterior"]][["a"]][["coeffs"]])

    for (draw in c(1, draws)) {
      expect_equal(fitted_of(back, draw), fitted_of(fitted, draw), info = paste("tvp:", tvp))
    }
    expect_equal(as.numeric(back[["data"]][["train"]][["w"]]),
                 as.numeric(centred_vec_model(tvp)[["data"]][["train"]][["w"]]))
    expect_null(attr(back[["data"]][["train"]][["w"]], "centre"))
    # The loadings and the cointegration vectors are not touched by centring.
    expect_equal(back[["posterior"]][["beta"]][["coeffs"]], fitted[["posterior"]][["beta"]][["coeffs"]])
    expect_error(rescale_error_correction(back), "attribute 'scale' or 'centre'")
  }
})

test_that("rescaling a centred model leaves its log-likelihood unchanged", {
  for (tvp in c(FALSE, TRUE)) {
    fitted <- centred_vec_fitted(tvp = tvp)
    back <- add_posterior_loglik(rescale_error_correction(fitted))

    expect_equal(as.numeric(back[["posterior"]][["loglik"]]),
                 as.numeric(fitted[["posterior"]][["loglik"]]),
                 info = paste("tvp:", tvp))
  }
})

test_that("rescaling a centred and scaled model keeps the fitted values", {
  fitted <- centred_vec_fitted(scale = TRUE)
  back <- rescale_error_correction(fitted)

  expect_equal(fitted_of(back, 1), fitted_of(fitted, 1))
  expect_equal(as.numeric(back[["data"]][["train"]][["w"]]),
               as.numeric(centred_vec_model()[["data"]][["train"]][["w"]]))
  expect_null(attr(back[["data"]][["train"]][["w"]], "scale"))
  expect_null(attr(back[["data"]][["train"]][["w"]], "centre"))
})

test_that("the starting values of the constant follow the centring", {
  for (tvp in c(FALSE, TRUE)) {
    model <- add_initial_values(centred_vec_priors(centred_vec_model(tvp), tvp))
    centred <- scale_error_correction(model, scale = FALSE, centre = TRUE)

    # The fit of the starting values, written as a draw of one row.
    as_draw <- function(object) {
      object[["posterior"]] <- list(
        a = list(coeffs = matrix(object[["initial"]][["a"]], 1)),
        beta = list(coeffs = matrix(object[["initial"]][["beta"]], 1)))
      object
    }

    expect_equal(fitted_of(as_draw(centred)), fitted_of(as_draw(model)), info = paste("tvp:", tvp))
    if (tvp) {
      expect_false(isTRUE(all.equal(centred[["initial"]][["a_init"]], model[["initial"]][["a_init"]])))
    }

    back <- rescale_error_correction(centred)
    expect_equal(back[["initial"]][["a"]], model[["initial"]][["a"]], info = paste("tvp:", tvp))
    expect_equal(back[["initial"]][["a_init"]], model[["initial"]][["a_init"]], info = paste("tvp:", tvp))
  }
})

test_that("vec_to_var and add_predictive_loglik refuse a centred model", {
  fitted <- centred_vec_fitted()

  expect_error(vec_to_var(fitted), "rescale_error_correction")

  # vec_to_var() recovers the levels from the error correction term, which needs
  # the means back.
  var <- vec_to_var(rescale_error_correction(fitted))
  expect_equal(as.numeric(var[["data"]][["train"]][["y"]]),
               as.numeric(stats::window(vec_data(), start = stats::start(var[["data"]][["train"]][["y"]]))))
})

test_that("bvec reconstructs the levels from a centred error correction term", {
  # bvec() takes the endogenous levels in 'w' and the deterministic regressors
  # in 'x_d', so the model has one lag and no restricted trend, as for the
  # scaled term.
  model <- create_bvecmodel(vec_data(), p = 1, r = 1, const = "unrestricted",
                            iterations = fx_iterations, burnin = fx_burnin)
  model <- scale_error_correction(model, centre = TRUE)
  set.seed(20260914)
  model <- add_initial_values(centred_vec_priors(model))
  fitted <- add_posterior_coefficients(model)
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

  expect_equal(as.numeric(object[["data"]][["original"]][["endogen"]]),
               as.numeric(stats::window(vec_data(),
                                        start = stats::start(object[["data"]][["train"]][["y"]]))))
})

test_that("the means survive a write and read round trip", {
  skip_if_not_installed("hdf5r")

  model <- scale_error_correction(centred_vec_model(), centre = TRUE)
  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(attr(restored[["data"]][["train"]][["w"]], "centre"),
               attr(model[["data"]][["train"]][["w"]], "centre"))
  expect_equal(unclass(rescale_error_correction(restored)[["data"]][["train"]][["w"]]),
               unclass(rescale_error_correction(model)[["data"]][["train"]][["w"]]),
               ignore_attr = TRUE)
})
