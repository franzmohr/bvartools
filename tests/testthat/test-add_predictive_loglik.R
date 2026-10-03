# One-step-ahead predictive densities of VAR and VEC models over expanding windows.

predlik_sigma <- function(error) {
  switch(error,
         wishart = list(df = "k", scale = 1),
         gamma = ,
         "gamma+covar" = list(shape = 3, rate = 0.01),
         list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
              state_variance = 0.05, offset = 1e-8))
}

predlik_vec <- function(tvp, error, r = 1, varsel = "none", iterations = 40) {
  model <- create_bvecmodel(vec_data(), p = 2, r = r, const = "unrestricted", tvp = tvp,
                            error = error, varsel = varsel,
                            iterations = iterations, burnin = 20)
  coef <- if (tvp) {
    list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4)
  } else {
    list(v_i = 1, v_i_det = 0.1)
  }
  coint <- if (tvp) {
    list(rho = 0.99, rho_min = 0.9, rho_max = 0.999)
  } else {
    list(v_i = 0, p_tau_i = 1)
  }
  args <- list(model, coef = coef, coint = coint, sigma = predlik_sigma(error))
  if (varsel == "bvs") {
    args[["varsel"]] <- list(inprior = 0.5, exclude_det = TRUE)
  }
  do.call(add_priors, args)
}

last_observation <- function(model) {
  .predictive_observation(model, .predictive_rank(model))
}

test_that("without innovations the density is the log-likelihood of the last period", {
  specs <- list(
    list(tvp = FALSE, error = "wishart"),
    list(tvp = FALSE, error = "gamma+covar"),
    list(tvp = FALSE, error = "sv+covar"),
    list(tvp = TRUE, error = "wishart"),
    list(tvp = TRUE, error = "gamma+covar"),
    list(tvp = TRUE, error = "sv"),
    list(tvp = TRUE, error = "sv+covar"),
    list(tvp = TRUE, error = "sv+covar", varsel = "bvs"),
    list(tvp = TRUE, error = "sv+covar", r = 0))

  for (spec in specs) {
    label <- paste(unlist(spec), collapse = " ")
    model <- do.call(predlik_vec, spec)
    set.seed(20110816)
    model <- add_initial_values(model)
    model <- add_posterior_loglik(add_posterior_coefficients(model))

    tt <- nrow(model[["data"]][["train"]][["y"]])
    density <- .predictive_log_density(model, last_observation(model), innovations = FALSE)
    expect_equal(density, as.numeric(model[["posterior"]][["loglik"]][, tt]),
                 tolerance = 1e-10, info = label)
  }
})

test_that("each window predicts the observation the next window adds", {
  model <- predlik_vec(tvp = FALSE, error = "wishart")
  windows <- use_expanding_window(model, start = stats::time(model[["data"]][["train"]][["y"]])[
    nrow(model[["data"]][["train"]][["y"]]) - 2])
  set.seed(1)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)
  windows <- add_predictive_loglik(windows)

  expect_s3_class(windows, "expandingwindow")
  n <- length(windows)
  expect_null(windows[[n]][["predictive"]])
  for (i in seq_len(n - 1)) {
    next_y <- windows[[i + 1]][["data"]][["train"]][["y"]]
    expect_equal(windows[[i]][["predictive"]][["period"]],
                 stats::time(next_y)[nrow(next_y)])
    expect_length(windows[[i]][["predictive"]][["loglik"]], 40L)
    # A constant model carries nothing forward, so the draws are deterministic.
    expect_equal(windows[[i]][["predictive"]][["loglik"]],
                 .predictive_log_density(windows[[i]], last_observation(windows[[i + 1]])))
  }
})

test_that("time varying states are carried forward before the density is taken", {
  set.seed(2)
  model <- add_initial_values(predlik_vec(tvp = TRUE, error = "sv+covar"))
  model <- add_posterior_coefficients(model)
  newdata <- last_observation(model)

  set.seed(3)
  first <- .predictive_log_density(model, newdata)
  set.seed(4)
  second <- .predictive_log_density(model, newdata)
  fixed <- .predictive_log_density(model, newdata, innovations = FALSE)

  expect_true(all(is.finite(first)))
  expect_false(isTRUE(all.equal(first, second)))
  expect_false(isTRUE(all.equal(first, fixed)))
})

test_that("the state variance of the log volatilities is read when it is stored", {
  set.seed(5)
  model <- add_initial_values(predlik_vec(tvp = TRUE, error = "sv"))
  model <- add_posterior_coefficients(model)
  k <- ncol(model[["data"]][["train"]][["y"]])
  tt <- nrow(model[["data"]][["train"]][["y"]])
  nd <- nrow(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  drawn <- .sv_state_variance(model, nd, k, tt)
  expect_identical(dim(drawn), c(nd, k))
  expect_true(all(drawn > 0))

  stored <- matrix(0.123, nd, k)
  model[["posterior"]][["u_sigma_inv"]][["sigma"]] <- stored
  expect_identical(.sv_state_variance(model, nd, k, tt), stored)
})

test_that("selection_criteria sums the predictive densities to the LPL", {
  models <- lapply(0:1, function(r) predlik_vec(tvp = FALSE, error = "wishart", r = r))
  class(models) <- c("modellist", "list")
  y <- models[[1]][["data"]][["train"]][["y"]]
  windows <- use_expanding_window(models, start = stats::time(y)[nrow(y) - 3])
  set.seed(6)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)
  windows <- add_predictive_loglik(windows)

  criteria <- selection_criteria(windows)
  expect_s3_class(criteria, "selcritlist")
  for (j in 1:2) {
    terms <- vapply(windows[[j]][-length(windows[[j]])], function(x) {
      .log_mean_exp(x[["predictive"]][["loglik"]])
    }, numeric(1))
    expect_equal(criteria[[j]][["LPL"]][["mean"]], sum(terms))
    expect_equal(attr(criteria[[j]][["LPL"]], "terms")[["lpd"]], unname(terms))
  }

  lpl <- vapply(criteria, function(x) x[["LPL"]][["mean"]], numeric(1))
  expect_identical(choose_best_model(criteria, criterion = "LPL"), which.max(lpl))
  expect_output(print(criteria), "Log predictive likelihood")
  expect_output(print(criteria[[1]]), "LPL")
})

test_that("windows the density cannot be taken for are refused", {
  # A VAR model is taken, so what a VAR window without draws is refused for is
  # the draws, as a VEC window would be.
  var_windows <- use_expanding_window(fx_var_model(), start = c(1997, 2))
  expect_error(add_predictive_loglik(var_windows), "no posterior draws")

  model <- predlik_vec(tvp = FALSE, error = "wishart")
  y <- model[["data"]][["train"]][["y"]]
  set.seed(7)
  windows <- add_initial_values(use_expanding_window(model, start = stats::time(y)[nrow(y) - 1]))
  expect_error(add_predictive_loglik(windows), "no posterior draws")

  windows <- add_posterior_coefficients(windows)
  # `[` drops the class of a list, so it is put back.
  single <- windows[1]
  class(single) <- class(windows)
  expect_error(add_predictive_loglik(single), "at least two windows")

  scaled <- windows
  attr(scaled[[1]][["data"]][["train"]][["w"]], "scale") <- rep(1, ncol(scaled[[1]][["data"]][["train"]][["w"]]))
  expect_error(add_predictive_loglik(scaled), "rescale_error_correction")

  gapped <- windows[c(1, 3)]
  class(gapped) <- class(windows)
  expect_error(add_predictive_loglik(gapped), "exactly one")
})

# --- VAR models -----------------------------------------------------------------
#
# A VAR model is the case of a VEC model without an error correction term, so the
# same two identities have to hold for it: without innovations the density is the
# pointwise log-likelihood of the last period, and every window predicts the
# observation the next window adds.

predlik_var <- function(tvp, error, varsel = "none", iterations = 40) {
  model <- create_bvarmodel(var_data(), p = 2, deterministic = "const", tvp = tvp,
                            error = error, varsel = varsel,
                            iterations = iterations, burnin = 20)
  coef <- if (tvp) {
    list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4)
  } else {
    list(v_i = 1, v_i_det = 0.1)
  }
  args <- list(model, coef = coef, sigma = predlik_sigma(error))
  if (varsel == "bvs") {
    args[["varsel"]] <- list(inprior = 0.5, exclude_det = TRUE)
  }
  do.call(add_priors, args)
}

test_that("a VAR model has a rank of zero in the density", {
  model <- predlik_var(tvp = FALSE, error = "wishart")
  expect_identical(.predictive_rank(model), 0L)
  expect_identical(.predictive_rank(predlik_vec(tvp = FALSE, error = "wishart")), 1L)
})

test_that("without innovations the density of a VAR is the log-likelihood of the last period", {
  specs <- list(
    list(tvp = FALSE, error = "wishart"),
    list(tvp = FALSE, error = "sv+covar"),
    list(tvp = TRUE, error = "wishart"),
    list(tvp = TRUE, error = "sv+covar"),
    list(tvp = TRUE, error = "sv+covar", varsel = "bvs"))

  for (spec in specs) {
    label <- paste(unlist(spec), collapse = " ")
    model <- do.call(predlik_var, spec)
    set.seed(20110816)
    model <- add_initial_values(model)
    model <- add_posterior_loglik(add_posterior_coefficients(model))

    tt <- nrow(model[["data"]][["train"]][["y"]])
    density <- .predictive_log_density(model, last_observation(model), innovations = FALSE)
    expect_equal(density, as.numeric(model[["posterior"]][["loglik"]][, tt]),
                 tolerance = 1e-10, info = label)
  }
})

test_that("the windows of a VAR model predict the observation the next one adds", {
  model <- predlik_var(tvp = FALSE, error = "wishart")
  y <- model[["data"]][["train"]][["y"]]
  windows <- use_expanding_window(model, start = stats::time(y)[nrow(y) - 2])
  set.seed(7)
  windows <- add_initial_values(windows)
  windows <- add_posterior_coefficients(windows)
  windows <- add_predictive_loglik(windows)

  n <- length(windows)
  expect_null(windows[[n]][["predictive"]])
  for (i in seq_len(n - 1)) {
    next_y <- windows[[i + 1]][["data"]][["train"]][["y"]]
    expect_equal(windows[[i]][["predictive"]][["period"]],
                 stats::time(next_y)[nrow(next_y)])
    expect_length(windows[[i]][["predictive"]][["loglik"]], 40L)
    expect_equal(windows[[i]][["predictive"]][["loglik"]],
                 .predictive_log_density(windows[[i]], last_observation(windows[[i + 1]])))
  }
  expect_equal(selection_criteria(windows)[["LPL"]][["mean"]],
               sum(vapply(windows[-n], function(x) {
                 .log_mean_exp(x[["predictive"]][["loglik"]])
               }, numeric(1))))
})

test_that("windows of mixed forms are refused", {
  vec <- predlik_vec(tvp = FALSE, error = "wishart")
  y <- vec[["data"]][["train"]][["y"]]
  windows <- use_expanding_window(vec, start = stats::time(y)[nrow(y) - 1])
  set.seed(8)
  windows <- add_posterior_coefficients(add_initial_values(windows))
  windows[[2]] <- structure(windows[[2]], class = c("list"))
  expect_error(add_predictive_loglik(windows), "VAR and VEC models only")
})

test_that("an algorithm whose likelihood is not normal is refused", {
  # A quantile VAR does store 'u_sigma_inv': the precision of the normal that
  # its scale mixture conditions on, period by period. The density would run on
  # it and return a number that is not the predictive density of the model, so
  # the algorithm is refused -- before the draws are looked at, which is why
  # these windows are never estimated.
  ald_windows <- function(tvp) {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                              error = "ald", quantile = 0.25, tvp = tvp,
                              iterations = fx_iterations, burnin = fx_burnin)
    y <- model[["data"]][["train"]][["y"]]
    use_expanding_window(model, start = stats::time(y)[nrow(y) - 1])
  }

  expect_error(add_predictive_loglik(ald_windows(tvp = FALSE)),
               "'VarNormalAld' of window 1")
  expect_error(add_predictive_loglik(ald_windows(tvp = TRUE)),
               "'VarTvpAld' of window 1")

  # A window that does not say what estimated it is refused as well: the
  # algorithms the density is the density of are named, and nothing else passes.
  unnamed <- ald_windows(tvp = FALSE)
  unnamed[[1]][["model"]][["algorithm"]] <- "VarSomethingNew"
  expect_error(add_predictive_loglik(unnamed), "'VarSomethingNew' of window 1")
  unnamed[[1]][["model"]][["algorithm"]] <- NULL
  expect_error(add_predictive_loglik(unnamed), "the algorithm of window 1")
})

test_that("a dummy that first enters during the evaluation period is left out of its forecast", {
  # use_expanding_window() leaves a dummy out of the windows that end before its
  # event, so the window just before the event has fewer regressors than the
  # observation the next window adds. Its density uses only its own regressors,
  # matched by name: the forecaster did not know about the event.
  n <- 80
  event <- 70
  freq <- 12
  set.seed(1)
  y <- stats::ts(matrix(cumsum(stats::rnorm(n)), ncol = 1, dimnames = list(NULL, "y")),
                 start = c(2018, 1), frequency = freq)
  y2 <- stats::ts(cbind(y1 = cumsum(stats::rnorm(n)), y2 = cumsum(stats::rnorm(n))),
                  start = c(2018, 1), frequency = freq)
  d <- stats::ts(matrix(as.numeric(seq_len(n + 24) == event), ncol = 1,
                        dimnames = list(NULL, "imp")),
                 start = c(2018, 1), frequency = freq)
  event_period <- stats::time(y)[event]

  specs <- list(
    list(label = "VAR, data dummy", vec = FALSE, dummy = list(data = d)),
    list(label = "VAR, impulse", vec = FALSE, dummy = list(impulse = c(2023, 10))),
    list(label = "VAR, step", vec = FALSE, dummy = list(step = c(2023, 10))),
    list(label = "VAR, p = 3, seasonal", vec = FALSE, p = 3, seasonal = TRUE,
         dummy = list(data = d)),
    list(label = "VEC, data dummy", vec = TRUE, dummy = list(data = d)))

  for (spec in specs) {
    p <- if (is.null(spec[["p"]])) 1 else spec[["p"]]
    model <- if (spec[["vec"]]) {
      create_bvecmodel(y2, p = p, r = 1, const = "unrestricted",
                       iterations = 30, burnin = 10)
    } else {
      create_bvarmodel(y, p = p, deterministic = "const",
                       seasonal = isTRUE(spec[["seasonal"]]),
                       iterations = 30, burnin = 10)
    }
    model <- do.call(add_dummy_variables, c(list(model), spec[["dummy"]]))
    model <- use_expanding_window(model, start = c(2023, 6))
    priors <- list(model, coef = list(v_i = 0.1, v_i_det = 0.01),
                   sigma = list(df = "k", scale = 1))
    if (spec[["vec"]]) {
      priors[["coint"]] <- list(v_i = 0, p_tau_i = 1)
    }
    model <- add_initial_values(do.call(add_priors, priors))
    model <- add_posterior_coefficients(model)
    model <- add_predictive_loglik(model)

    n_windows <- length(model)
    expect_null(model[[n_windows]][["predictive"]], info = spec[["label"]])
    periods <- numeric(n_windows - 1)
    for (i in seq_len(n_windows - 1)) {
      expect_false(is.null(model[[i]][["predictive"]]), info = spec[["label"]])
      periods[i] <- model[[i]][["predictive"]][["period"]]
    }

    i <- which(abs(periods - event_period) < 1e-8)
    expect_length(i, 1)
    # The window that predicts the event does not have the dummy; the next does.
    dummy_name <- model[[i + 1]][["model"]][["dummy_variables"]][["name"]]
    expect_false(any(dummy_name %in% colnames(model[[i]][["data"]][["train"]][["x"]])),
                 info = spec[["label"]])
    expect_true(all(dummy_name %in% colnames(model[[i + 1]][["data"]][["train"]][["x"]])),
                info = spec[["label"]])
    loglik <- model[[i]][["predictive"]][["loglik"]]
    expect_true(all(is.finite(loglik)), info = spec[["label"]])

    # It is the density of the observation without the dummy.
    following <- model[[i + 1]]
    x_names <- colnames(model[[i]][["data"]][["train"]][["x"]])
    x_next <- following[["data"]][["train"]][["x"]]
    newdata <- last_observation(following)
    newdata[["x"]] <- as.numeric(x_next[nrow(x_next), x_names])
    expect_equal(loglik, .predictive_log_density(model[[i]], newdata),
                 info = spec[["label"]])
  }
})
