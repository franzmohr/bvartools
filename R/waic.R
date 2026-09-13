# Widely applicable information criterion.
#
# AIC, BIC and HQ penalise a model by the number of parameters it nominally
# has. That number says little about a model whose coefficients or variances
# follow a state equation: a time varying parameter model carries a state per
# coefficient and period, but the prior on the state variances decides how much
# of that freedom is actually used, so the count of parameters and the
# flexibility of the fit part company. WAIC penalises instead by the variance
# of the pointwise log-likelihood across draws, which measures the flexibility
# that was used, and is therefore the criterion to reach for when constant,
# time varying and stochastic volatility specifications are compared with each
# other.
#
# Following Watanabe (2010) in the form given by Vehtari, Gelman and Gabry
# (2017), with the pointwise log-likelihood of a draws x periods matrix,
#
#   lppd   = sum_t log( mean_s exp( loglik[s, t] ) )
#   p_waic = sum_t var_s( loglik[s, t] )
#   WAIC   = -2 (lppd - p_waic),
#
# so that WAIC is on the deviance scale of the other criteria and smaller is
# better. The standard error is the one of a sum of T independent terms, which
# is what makes a difference between two models readable: a gap smaller than a
# few standard errors is not one to choose on.
.waic <- function(loglik) {

  loglik <- as.matrix(loglik)
  draws <- nrow(loglik)
  tt <- ncol(loglik)

  if (draws < 2) {
    return(NULL)
  }

  # log(mean(exp(z))) with the maximum factored out, since exp() of a
  # log-likelihood in the hundreds overflows long before the mean is taken.
  log_mean_exp <- function(z) {
    z_max <- max(z)
    z_max + log(mean(exp(z - z_max)))
  }

  lppd <- apply(loglik, 2, log_mean_exp)
  p_waic <- apply(loglik, 2, stats::var)

  elpd <- lppd - p_waic
  waic <- -2 * sum(elpd)

  # On the deviance scale, hence the factor of two.
  se <- 2 * sqrt(tt * stats::var(elpd))

  # Where the variance of the pointwise log-likelihood of a period exceeds 0.4,
  # the correction it contributes cannot be relied on (Vehtari, Gelman and Gabry,
  # 2017). The count of such periods travels with the criterion.
  var_threshold <- 0.4

  list("waic" = waic, "se" = se, "p_waic" = sum(p_waic),
       "n_high_var" = sum(p_waic > var_threshold),
       "var_threshold" = var_threshold)
}

# The criteria are reported as data frames of the same four columns, and
# choose_best_model() reads the 'mean' one. WAIC is a point estimate rather
# than a quantity with a posterior distribution, so 'mean' and 'median' hold
# the estimate and the two quantile columns hold a normal interval built from
# its standard error, which is the usual way it is reported.
.waic_data_frame <- function(loglik, ci_low, ci_high) {

  result <- .waic(loglik)
  if (is.null(result)) {
    return(NULL)
  }

  z <- stats::qnorm(ci_high)

  out <- data.frame("mean" = result[["waic"]],
                    "median" = result[["waic"]],
                    "qlower" = result[["waic"]] - z * result[["se"]],
                    "qupper" = result[["waic"]] + z * result[["se"]])

  # The diagnostics travel as attributes, as those of LOOIC do.
  attr(out, "p_waic") <- result[["p_waic"]]
  attr(out, "n_high_var") <- result[["n_high_var"]]
  attr(out, "var_threshold") <- result[["var_threshold"]]

  out
}

# A note under the table of criteria when the variance of the pointwise
# log-likelihood is too large in some periods for the correction built from it
# to be relied on. That variance is WAIC's penalty, so a time varying or
# stochastic volatility model sampled briefly can report a WAIC far below its
# fit with nothing in the number itself to say so. AIC, BIC and HQ start from
# the deviance at the posterior mean instead and are not affected.
.print_waic_diagnostics <- function(waic) {

  if (is.null(waic)) {
    return(invisible(NULL))
  }

  n_high_var <- attr(waic, "n_high_var")
  if (is.null(n_high_var) || n_high_var == 0) {
    return(invisible(NULL))
  }

  cat("\nWAIC: ", n_high_var,
      ifelse(n_high_var == 1, " period has", " periods have"),
      " a pointwise log-likelihood variance above ",
      format(attr(waic, "var_threshold")),
      ", so the correction it rests on is unreliable there.\n", sep = "")

  invisible(NULL)
}

# The same note for a list of models, naming the models concerned rather than
# repeating a line per model.
.print_waic_diagnostics_list <- function(x) {

  n_high_var <- vapply(x, function(y) {
    value <- attr(y[["WAIC"]], "n_high_var")
    if (is.null(value)) 0L else as.integer(value)
  }, integer(1))

  if (all(n_high_var == 0)) {
    return(invisible(NULL))
  }

  affected <- which(n_high_var > 0)
  threshold <- attr(x[[affected[1]]][["WAIC"]], "var_threshold")

  cat("\nPeriods with a pointwise log-likelihood variance above ", format(threshold),
      ", which makes the correction of WAIC unreliable, in model",
      ifelse(length(affected) == 1, " ", "s "),
      paste0(affected, " (", n_high_var[affected], ")", collapse = ", "), ".\n", sep = "")

  invisible(NULL)
}
