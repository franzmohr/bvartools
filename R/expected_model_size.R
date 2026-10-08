
# The size of every block of draws a VAR or VEC model will hold, from its
# specification alone. The numbers are the allocations of the samplers in
# src/core/models/ and of their bindings in src/*.cpp: every block is a matrix
# of doubles with one row per kept draw, and what decides its width is written
# beside it. test-expected_model_size.R holds each of them against a fitted model, so
# a sampler that stores something new fails there rather than here.
.model_size <- function(object, chains = NULL) {

  model <- object[["model"]]
  priors <- object[["priors"]]
  k <- model[["k"]]
  tt <- NROW(object[["data"]][["train"]][["y"]])
  chains <- .check_chains(object, chains)

  elements <- character(0)
  steps <- character(0)
  draws <- numeric(0)
  columns <- numeric(0)
  add <- function(element, step, n_draws, n_columns) {
    if (n_draws > 0 && n_columns > 0) {
      elements <<- c(elements, element)
      steps <<- c(steps, step)
      draws <<- c(draws, n_draws)
      columns <<- c(columns, n_columns)
    }
  }
  coef_step <- "add_posterior_coefficients"
  vec <- inherits(object, "bvecmodel")
  rank <- if (vec) model[["rank"]] else 0
  n_beta <- if (vec) model[["k_beta"]] * rank else 0

  if (.is_discount(object)) {

    # Not a chain: one row per period of a closed form, and the fixed
    # cointegration matrix of a VEC as a single row.
    n_design <- NCOL(object[["data"]][["train"]][["x"]]) + rank
    add("posterior$a$mean", coef_step, tt, n_design * k)
    add("posterior$a$scale", coef_step, tt, n_design * k)
    add("posterior$a$cov", coef_step, tt, n_design^2)
    add("posterior$u_sigma$scale", coef_step, tt, k^2)
    add("posterior$df", coef_step, tt, 1)
    add("posterior$beta$coeffs", coef_step, 1, n_beta)
    add("posterior$loglik", "add_posterior_loglik", 1, tt)

  } else {

    n <- as.numeric(model[["iterations"]]) * chains
    tvp <- isTRUE(model[["tvp"]])
    periods <- if (tvp) tt else 1
    error <- model[["error"]]
    if (is.null(error)) {
      error <- "wishart"
    }
    sv <- grepl("^sv", error)
    ald <- error == "ald"
    covar <- grepl("covar", error) && k > 1

    # The draws of the non-centred parameterisation, beside the state variances
    # of a block whose prior is 'omega_v'. The prior of the log-volatilities is
    # 'u_sigma' and their draws are kept under 'u_sigma_inv'.
    noncentred <- function(block, width, element = block) {
      if (tvp && !is.null(priors[[block]][["omega_v"]])) {
        add(paste0("posterior$", element, "$omega"), coef_step, n, width)
        add(paste0("posterior$", element, "$omega_log_zero"), coef_step, n, width)
        add(paste0("posterior$", element, "$omega_log_zero_joint"), coef_step, n, 1)
      }
    }

    # The coefficients: one per column of the SUR design, which already holds
    # the loadings of a VEC and the contemporaneous terms of a structural model.
    n_a <- NCOL(object[["data"]][["train"]][["z"]])
    add("posterior$a$coeffs", coef_step, n, n_a * periods)
    if (tvp) {
      add("posterior$a$sigma", coef_step, n, n_a)
    }
    if (!is.null(model[["varsel"]]) && model[["varsel"]] != "none") {
      add("posterior$a$lambda", coef_step, n, n_a)
    }
    if (n_a > 0) {
      noncentred("a", n_a)
    }

    add("posterior$beta$coeffs", coef_step, n, n_beta * periods)
    if (tvp && !is.null(priors[["beta"]][["rho_min"]])) {
      add("posterior$beta$rho", coef_step, n, 1)
    }

    if (covar) {
      n_psi <- k * (k - 1) / 2
      add("posterior$psi$coeffs", coef_step, n, k^2 * periods)
      if (tvp) {
        add("posterior$psi$sigma", coef_step, n, n_psi)
      }
      if (!is.null(priors[["psi"]][["varsel"]]) && priors[["psi"]][["varsel"]] != "none") {
        add("posterior$psi$lambda", coef_step, n, k^2)
      }
      noncentred("psi", n_psi)
    }

    # The error precision moves with time under stochastic volatility and the
    # asymmetric Laplace, and with a covariance block that drifts.
    add("posterior$u_sigma_inv$coeffs", coef_step, n,
        k^2 * (if (sv || ald || (tvp && covar)) tt else 1))
    if (sv) {
      add("posterior$u_sigma_inv$sigma", coef_step, n, k)
      noncentred("u_sigma", k, "u_sigma_inv")
    }
    if (error != "wishart") {
      add("posterior$u_omega_inv$coeffs", coef_step, n, k * (if (sv || ald) tt else 1))
    }
    if (ald) {
      add("posterior$u_scale$coeffs", coef_step, n, k)
    }

    add("posterior$loglik", "add_posterior_loglik", n, tt)
  }

  # Forecast draws are counted once add_forecast_input() has set a horizon.
  if (!is.null(model[["h"]])) {
    n_forecast <- as.numeric(model[["iterations"]]) * chains
    add("posterior$forecast$forecasts", "add_posterior_forecasts", n_forecast, k * model[["h"]])
  }

  result <- data.frame(element = c("data, priors and initial values", elements),
                       step = c("", steps),
                       draws = c(NA, draws),
                       columns = c(NA, columns),
                       bytes = c(as.numeric(utils::object.size(object[names(object) != "posterior"])),
                                 8 * draws * columns),
                       stringsAsFactors = FALSE)
  class(result) <- c("modelsize", "data.frame")
  result
}
