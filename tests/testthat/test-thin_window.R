test_that("thin keeps every nth draw", {
  model <- fx_var_fitted()
  thinned <- thin(model, thin = 3)
  draws <- model[["posterior"]][["a"]][["coeffs"]]

  expect_s3_class(thinned, "bvarmodel")
  expect_identical(nrow(thinned[["posterior"]][["a"]][["coeffs"]]),
                   as.integer(fx_iterations / 3))
  expect_equal(unname(as.matrix(thinned[["posterior"]][["a"]][["coeffs"]])),
               unname(as.matrix(draws[seq(3, fx_iterations, 3), ])))
})

test_that("thin applies to every posterior element", {
  thinned <- thin(fx_var_fitted(), thin = 2)
  expected <- as.integer(fx_iterations / 2)

  expect_identical(nrow(thinned[["posterior"]][["a"]][["coeffs"]]), expected)
  expect_identical(nrow(thinned[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
                   expected)
  expect_identical(nrow(thinned[["posterior"]][["loglik"]]), expected)
})

test_that("thinning by one leaves the draws untouched", {
  expect_equal(thin(fx_var_fitted(), thin = 1)[["posterior"]],
               fx_var_fitted()[["posterior"]])
})

test_that("thin works for VEC models and model lists", {
  vec <- thin(fx_vec_fitted(), thin = 2)
  models <- thin(fx_var_modellist(), thin = 2)

  expect_identical(nrow(vec[["posterior"]][["beta"]][["coeffs"]]),
                   as.integer(fx_iterations / 2))
  expect_s3_class(models, "modellist")
  expect_true(all(vapply(models,
                         function(x) nrow(x[["posterior"]][["a"]][["coeffs"]]),
                         integer(1)) == fx_iterations / 2))
})

test_that("window subsets the estimation sample in time", {
  model <- stats::window(fx_var_model(), start = c(1984, 2), end = c(1990, 1))
  y <- model[["data"]][["train"]][["y"]]

  expect_equal(stats::tsp(y)[1], 1984.25)
  expect_equal(stats::tsp(y)[2], 1990)
  expect_identical(nrow(y), 24L)
  # x is cut alongside y.
  expect_identical(nrow(model[["data"]][["train"]][["x"]]), nrow(y))
})

test_that("window cuts the SUR matrix in matching row blocks", {
  full <- fx_var_model()
  cut <- stats::window(full, start = c(1984, 2), end = c(1990, 1))
  k <- full[["model"]][["k"]]

  expect_identical(nrow(cut[["data"]][["train"]][["z"]]),
                   as.integer(k * nrow(cut[["data"]][["train"]][["y"]])))

  # The retained rows are the ones belonging to the retained periods.
  keep <- which(stats::time(full[["data"]][["train"]][["y"]]) %in%
                  stats::time(cut[["data"]][["train"]][["y"]]))
  rows <- rep((keep - 1) * k, each = k) + seq_len(k)
  expect_equal(cut[["data"]][["train"]][["z"]],
               full[["data"]][["train"]][["z"]][rows, ])
})

test_that("window works for VEC models and model lists", {
  vec <- stats::window(fx_vec_model(), start = c(1987, 1))
  models <- stats::window(fx_var_modellist(), start = c(1984, 2))

  expect_equal(stats::tsp(vec[["data"]][["train"]][["y"]])[1], 1987)
  expect_s3_class(models, "modellist")
  expect_true(all(vapply(
    models,
    function(x) stats::tsp(x[["data"]][["train"]][["y"]])[1],
    numeric(1)
  ) == 1984.25))
})

test_that("thin keeps every block of draws in step", {
  keep <- seq(3, fx_iterations, 3)

  for (model in list(fx_var_tvp_fitted("gamma"), fx_vec_tvp_fitted("sv"))) {
    before <- draw_blocks(model[["posterior"]])
    after <- draw_blocks(thin(model, thin = 3)[["posterior"]])

    expect_named(after, names(before))
    for (name in names(before)) {
      expect_equal(unname(as.matrix(after[[name]])),
                   unname(as.matrix(before[[name]])[keep, , drop = FALSE]), label = name)
      expect_equal(attr(after[[name]], "mcpar"), c(3, max(keep), 3), label = name)
    }
  }

  # The blocks a list of names used to miss are among the ones checked.
  expect_true(all(c("a$sigma", "beta$rho") %in%
                    names(draw_blocks(fx_vec_tvp_fitted("sv")[["posterior"]]))))
})

test_that("window cuts the posterior paths to the periods it keeps", {
  specs <- list(list(model = fx_var_tvp_fitted("sv"), start = c(1984, 2), end = c(1995, 1)),
                list(model = fx_vec_tvp_fitted("sv"), start = c(1987, 1), end = c(1997, 4)))

  for (spec in specs) {
    full <- spec[["model"]]
    cut <- stats::window(full, start = spec[["start"]], end = spec[["end"]])
    k <- full[["model"]][["k"]]
    keep <- which(stats::time(full[["data"]][["train"]][["y"]]) %in%
                    stats::time(cut[["data"]][["train"]][["y"]]))
    columns <- function(width) rep((keep - 1) * width, each = width) + seq_len(width)

    widths <- c("a$coeffs" = ncol(full[["data"]][["train"]][["z"]]),
                "u_sigma_inv$coeffs" = k^2, "u_omega_inv$coeffs" = k, "loglik" = 1)
    if (!is.null(full[["posterior"]][["beta"]])) {
      widths[["beta$coeffs"]] <- ncol(full[["data"]][["train"]][["w"]]) * full[["model"]][["rank"]]
    }

    before <- draw_blocks(full[["posterior"]])
    after <- draw_blocks(cut[["posterior"]])
    for (name in names(widths)) {
      expect_equal(unname(as.matrix(after[[name]])),
                   unname(as.matrix(before[[name]])[, columns(widths[[name]]), drop = FALSE]),
                   label = name)
      expect_equal(attr(after[[name]], "mcpar"), attr(before[[name]], "mcpar"), label = name)
    }

    # Draws that do not vary by period are left as they are.
    for (name in setdiff(names(before), names(widths))) {
      expect_identical(after[[name]], before[[name]], label = name)
    }

    # A period of the window is the period of the original sample it came from.
    expect_equal(summary(cut)[["a"]][["means"]],
                 summary(full, period = max(keep))[["a"]][["means"]])
    expect_equal(summary(cut, period = 1)[["sigma"]][["means"]],
                 summary(full, period = min(keep))[["sigma"]][["means"]])
  }
})

test_that("the thinning interval is a positive integer within the draws", {
  model <- fx_at_vec_tvp()

  expect_error(thin(model, thin = fx_iterations + 1),
               paste0("only ", fx_iterations, " draws"))
  expect_error(thin(model, thin = 0), "single positive integer")
  expect_error(thin(model, thin = 2.5), "single positive integer")
  expect_error(thin(model, thin = c(2, 3)), "single positive integer")
  expect_identical(nrow(thin(model, thin = fx_iterations)[["posterior"]][["a"]][["coeffs"]]), 1L)
})

# A model whose sampler keeps one draw in `thin`, built from a fitted fixture's
# starting point so that its chain and the unthinned one start from the same place.
with_sampler_thin <- function(model, thin) {
  model[["model"]][["iterations"]] <- as.integer(fx_iterations / thin)
  model[["model"]][["thin"]] <- as.integer(thin)
  model
}

test_that("create_*model() takes the thinning interval of the sampler", {
  expect_null(fx_var_model()[["model"]][["thin"]])
  expect_null(create_bvarmodel(var_data(), p = 1, iterations = 10, burnin = 5,
                               thin = 1)[["model"]][["thin"]])
  expect_identical(create_bvarmodel(var_data(), p = 1, iterations = 10, burnin = 5,
                                    thin = 4)[["model"]][["thin"]], 4L)
  expect_identical(create_bvecmodel(vec_data(), p = 1, r = 1, const = "unrestricted",
                                    iterations = 10, burnin = 5, thin = 2)[["model"]][["thin"]], 2L)
  expect_error(create_bvarmodel(var_data(), p = 1, thin = 0), "single positive integer")
  expect_error(create_bvecmodel(vec_data(), p = 1, r = 1, thin = 2.5), "single positive integer")
})

test_that("a thinned chain is every thin-th draw of the unthinned one", {
  keep <- seq(3, fx_iterations, 3)

  for (base in list(fx_var_initial(), fx_vec_initial())) {
    set.seed(20260913)
    full <- add_posterior_coefficients(base)
    set.seed(20260913)
    thinned <- add_posterior_coefficients(with_sampler_thin(base, 3))

    before <- draw_blocks(full[["posterior"]])
    after <- draw_blocks(thinned[["posterior"]])
    expect_named(after, names(before))
    for (name in names(before)) {
      # Exactly: which draws are kept is the only thing thinning may change.
      expect_identical(c(as.matrix(after[[name]])),
                       c(as.matrix(before[[name]])[keep, , drop = FALSE]), label = name)
      expect_equal(attr(after[[name]], "mcpar"), c(3, fx_iterations, 3), label = name)
    }
  }

  expect_identical(vec_to_var(add_posterior_coefficients(
    with_sampler_thin(fx_vec_initial(), 3)))[["model"]][["thin"]], 3L)
})

test_that("the log likelihood and the forecasts of a thinned chain carry its labels", {
  set.seed(7)
  model <- add_posterior_coefficients(with_sampler_thin(fx_var_initial(), 3))
  model <- add_posterior_loglik(model)
  model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 2))

  for (block in c("a", "loglik", "forecast")) {
    draws <- switch(block,
                    a = model[["posterior"]][["a"]][["coeffs"]],
                    forecast = model[["posterior"]][["forecast"]][["forecasts"]],
                    model[["posterior"]][[block]])
    expect_equal(attr(draws, "mcpar"), c(3, fx_iterations, 3), label = block)
  }
})

test_that("thin() on a thinned chain counts from the chain's own labels", {
  set.seed(8)
  model <- add_posterior_coefficients(with_sampler_thin(fx_var_initial(), 3))
  thinned <- thin(model, thin = 2)
  draws <- thinned[["posterior"]][["a"]][["coeffs"]]

  expect_identical(nrow(draws), 5L)
  expect_equal(attr(draws, "mcpar"), c(6, fx_iterations, 6))
  expect_identical(c(as.matrix(draws)),
                   c(as.matrix(model[["posterior"]][["a"]][["coeffs"]])[seq(2, 10, 2), ]))
})

test_that("window() cuts the Psi path of a time varying covariance", {
  # The samplers store the whole k x k Psi per period. window() used to expect
  # its k(k - 1)/2 free elements, so the path was left at the full sample and a
  # forecast or score from the window read Psi from a period of the original.
  model <- fx_var_tvp_fitted("gamma+covar")
  k <- model[["model"]][["k"]]
  y <- model[["data"]][["train"]][["y"]]
  tt <- nrow(y)
  w <- window(model, start = stats::time(y)[11])

  expect_identical(ncol(w[["posterior"]][["psi"]][["coeffs"]]), as.integer(k * k * (tt - 10)))
  expect_identical(ncol(w[["posterior"]][["u_sigma_inv"]][["coeffs"]]), as.integer(k * k * (tt - 10)))
  # Cutting from the start leaves the last period where it was.
  last <- function(x, n) {
    x <- as.matrix(x)
    x[, ncol(x) - n + seq_len(n), drop = FALSE]
  }
  expect_equal(last(w[["posterior"]][["psi"]][["coeffs"]], k * k),
               last(model[["posterior"]][["psi"]][["coeffs"]], k * k), ignore_attr = TRUE)
})

test_that("window() cuts the posterior of a discounted model to its periods", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            algorithm = "discount", delta_beta = 0.99, delta_sigma = 0.98,
                            iterations = 10, burnin = 0)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      sigma = list(df = "k", scale = 1))
  model <- add_posterior_coefficients(add_initial_values(model))
  y <- model[["data"]][["train"]][["y"]]
  keep <- nrow(y) - 5L
  w <- window(model, end = stats::time(y)[keep])

  for (block in list(c("a", "mean"), c("a", "scale"), c("a", "cov"),
                     c("u_sigma", "scale"), "df")) {
    what <- paste(block, collapse = "$")
    expect_identical(NROW(w[["posterior"]][[block]]), keep, info = what)
    expect_equal(unname(as.matrix(w[["posterior"]][[block]])),
                 unname(as.matrix(model[["posterior"]][[block]])[seq_len(keep), , drop = FALSE]),
                 info = what)
  }
})

test_that("pooled chains are thinned one chain at a time", {
  # Thinned as one sequence, the second chain kept different positions than
  # the first whenever a chain's length was not a multiple of 'thin'.
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = 17, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
                      sigma = list(df = "k", scale = 1))
  model <- add_posterior_coefficients(add_initial_values(model), chains = 2)
  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])

  thinned <- thin(model, thin = 4)
  kept <- as.matrix(thinned[["posterior"]][["a"]][["coeffs"]])
  expect_identical(nrow(kept), 8L)
  expect_equal(unname(kept), unname(a[c(4, 8, 12, 16, 21, 25, 29, 33), ]))
  expect_no_error(chain_diagnostics(thinned))
})

test_that("a discounted model says there is nothing to thin", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            algorithm = "discount", iterations = 10, burnin = 0)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      sigma = list(df = "k", scale = 1))
  model <- add_posterior_coefficients(add_initial_values(model))
  expect_error(thin(model, thin = 2), "nothing to thin")
})
