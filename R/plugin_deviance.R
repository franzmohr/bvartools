# Deviance of a model at its point estimate.
#
# AIC, BIC and HQ are not penalties bolted onto an arbitrary measure of fit.
# Each estimates how the model would do on data it has not seen, and each
# corrects the deviance of a *fitted* model for the optimism of having fitted
# it. The correction is derived for a likelihood that was maximised, and it is
# the deviance at that maximum -- or, here, at the Bayesian point estimate --
# that the correction belongs to. The mean of the deviance over the posterior
# exceeds it by about the effective number of parameters, so a penalty added to
# the mean charges the complexity of the model twice.
#
# The deviance at the point estimate is evaluated directly: the model's own
# log-likelihood is computed once at the posterior mean of its parameters. It
# used to be recovered from the draws of the log-likelihood as their mean
# deviance less the sum of the pointwise variances, which is WAIC's penalty.
# That sum is not the effective number of parameters. On the Austrian data of
# at_macrodata, VARs with a flat prior and lag orders 0 to 4 have 9 to 45
# parameters, and the sum came to 11 to 58 -- a quarter to a third more -- while
# the deviance at the posterior mean agrees with the maximum likelihood one. So
# AIC, BIC and HQ came out too low by an amount that grew with the size of the
# model, which reordered lag orders under AIC.
#
# Evaluating the model at a point works for every specification the package
# estimates, time varying and stochastic volatility models included, because
# add_posterior_loglik() already knows how to score each of them. What has to be
# decided is only what the point is. For every block of draws it is the mean.
# The loadings and cointegration vectors of a VEC model are the exception:
# alpha and beta are identified only up to a rotation, so their means say
# nothing, and the mean of Pi = alpha beta' is generally of full rank. The point
# is therefore the best approximation of rank r to the posterior mean of Pi,
# factored back into alpha and beta.
.plugin_deviance <- function(object) {

  point <- .posterior_mean_model(object)

  loglik <- tryCatch(add_posterior_loglik(point)[["posterior"]][["loglik"]],
                     error = function(e) {
                       warning("The log-likelihood could not be evaluated at the posterior mean, ",
                               "so AIC, BIC and HQ are not available: ", conditionMessage(e),
                               call. = FALSE)
                       NULL
                     })

  if (is.null(loglik)) {
    return(NA_real_)
  }

  -2 * sum(loglik)
}

# The model with its posterior replaced by a single draw, the posterior mean.
.posterior_mean_model <- function(object) {

  posterior <- object[["posterior"]]

  # Blocks that are not parameters of the likelihood, or whose mean is not a
  # point of the model: forecasts, their errors, the log-likelihood itself and
  # the rotations of a sign restricted identification.
  skip <- c("loglik", "forecast", "forecast_errors", "q")

  one_draw <- function(draws) {
    draws <- as.matrix(draws)
    coda::mcmc(matrix(colMeans(draws), nrow = 1), start = 1, end = 1, thin = 1)
  }

  point <- list()
  for (name in setdiff(names(posterior), skip)) {
    block <- posterior[[name]]
    if (is.list(block)) {
      point[[name]] <- lapply(block, function(draws) {
        if (is.null(draws)) NULL else one_draw(draws)
      })
    }
  }

  rank <- object[["model"]][["rank"]]
  if (inherits(object, "bvecmodel") && isTRUE(rank > 0) &&
      !is.null(posterior[["beta"]][["coeffs"]])) {
    point <- .posterior_mean_pi(object, point)
  }

  object[["posterior"]] <- point
  object[["model"]][["iterations"]] <- 1L
  object[["model"]][["burnin"]] <- 0L
  object[["model"]][["thin"]] <- NULL

  object
}

# Replaces alpha and beta in the posterior mean of a VEC model by the rank r
# factorisation of the posterior mean of Pi, period by period for a time varying
# model. 'a' holds alpha in its first k r elements per period, as a k x r matrix
# stacked by column, and 'beta' holds the k_beta x r matrix of the period.
.posterior_mean_pi <- function(object, point) {

  k <- object[["model"]][["k"]]
  rank <- object[["model"]][["rank"]]
  k_beta <- ncol(object[["data"]][["train"]][["w"]])
  n_alpha <- k * rank
  n_beta <- k_beta * rank

  a <- as.matrix(object[["posterior"]][["a"]][["coeffs"]])
  beta <- as.matrix(object[["posterior"]][["beta"]][["coeffs"]])
  draws <- nrow(a)
  periods <- ncol(beta) / n_beta
  n_a <- ncol(a) / periods

  a_point <- as.numeric(point[["a"]][["coeffs"]])
  beta_point <- as.numeric(point[["beta"]][["coeffs"]])

  for (t in seq_len(periods)) {
    pos_a <- (t - 1) * n_a
    pos_beta <- (t - 1) * n_beta

    # E[alpha beta'] summed over the r columns, without forming a matrix per draw
    pi_mean <- matrix(0, k, k_beta)
    for (q in seq_len(rank)) {
      alpha_q <- a[, pos_a + (q - 1) * k + seq_len(k), drop = FALSE]
      beta_q <- beta[, pos_beta + (q - 1) * k_beta + seq_len(k_beta), drop = FALSE]
      pi_mean <- pi_mean + crossprod(alpha_q, beta_q) / draws
    }

    decomposition <- svd(pi_mean, nu = rank, nv = rank)
    alpha_t <- decomposition[["u"]] %*% diag(decomposition[["d"]][seq_len(rank)], rank)
    beta_t <- decomposition[["v"]]

    a_point[pos_a + seq_len(n_alpha)] <- as.numeric(alpha_t)
    beta_point[pos_beta + seq_len(n_beta)] <- as.numeric(beta_t)
  }

  point[["a"]][["coeffs"]] <- coda::mcmc(matrix(a_point, nrow = 1), start = 1, end = 1, thin = 1)
  point[["beta"]][["coeffs"]] <- coda::mcmc(matrix(beta_point, nrow = 1), start = 1, end = 1, thin = 1)

  point
}

# The criteria that are point estimates rather than quantities with a
# posterior distribution are reported in the same four columns as the others,
# with the estimate in both the mean and the median and no credible band,
# which they do not have: a criterion is a function of the data and the
# estimator, not a parameter to be uncertain about.
.point_criterion <- function(value) {
  data.frame("mean" = value,
             "median" = value,
             "qlower" = NA_real_,
             "qupper" = NA_real_)
}
