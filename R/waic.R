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

  list("waic" = waic, "se" = se, "p_waic" = sum(p_waic))
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

  data.frame("mean" = result[["waic"]],
             "median" = result[["waic"]],
             "qlower" = result[["waic"]] - z * result[["se"]],
             "qupper" = result[["waic"]] + z * result[["se"]])
}
