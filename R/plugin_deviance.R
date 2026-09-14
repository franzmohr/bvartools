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
# is therefore the approximation of rank r to the posterior mean of Pi that
# fits best in the metric of the likelihood, factored back into alpha and beta;
# see .posterior_mean_pi().
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

# Replaces alpha and beta in the posterior mean of a VEC model by a rank r
# factorisation of the posterior mean of Pi, period by period for a time varying
# model. 'a' holds alpha in its first k r elements per period, as a k x r matrix
# stacked by column, and 'beta' holds the k_beta x r matrix of the period.
#
# Which matrix of rank r is closest depends on how a difference in Pi is
# measured, and it has to be measured as the likelihood measures it: a change D
# in Pi moves the errors of period t by D w_t, which the likelihood weighs with
# the error precision. The point therefore minimises
#
#     sum_t (D w_t)' Q (D w_t) = tr(Q^(1/2) D (W'W) D' Q^(1/2)),  D = Pi_mean - Pi_r,
#
# with Q the posterior mean of the error precision (of the period, if it moves)
# and W the error correction term. Its solution is the truncated SVD of
# Q^(1/2) Pi_mean (W'W)^(1/2), transformed back. The plain truncated SVD of
# Pi_mean, which was used before, measures every element of Pi alike. The series
# in W are levels on scales of their own -- output times 100 is about 450 in
# at_macrodata, an interest rate about 2 -- so it gave up fit where the series
# are large to keep elements that hardly matter: in a four-variable VEC model of
# at_macrodata levels the AIC of rank 2 came to 2213 against -121 at the maximum
# likelihood estimate, and AIC ranked the ranks in an order unrelated to theirs.
# The weighted point is also invariant to rescaling the series or the equations.
# A time varying model is approximated period by period in the metric of the
# whole sample's W'W.
.posterior_mean_pi <- function(object, point) {

  k <- object[["model"]][["k"]]
  rank <- object[["model"]][["rank"]]
  w <- as.matrix(object[["data"]][["train"]][["w"]])
  k_beta <- ncol(w)
  tt <- nrow(w)
  n_alpha <- k * rank
  n_beta <- k_beta * rank

  a <- as.matrix(object[["posterior"]][["a"]][["coeffs"]])
  beta <- as.matrix(object[["posterior"]][["beta"]][["coeffs"]])
  draws <- nrow(a)
  periods <- ncol(beta) / n_beta
  n_a <- ncol(a) / periods

  a_point <- as.numeric(point[["a"]][["coeffs"]])
  beta_point <- as.numeric(point[["beta"]][["coeffs"]])

  metric_w <- .symmetric_roots(crossprod(w))
  u_sigma_inv <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]]
  precision_mean <- if (is.null(u_sigma_inv)) NULL else colMeans(as.matrix(u_sigma_inv))
  precision_path <- !is.null(precision_mean) && length(precision_mean) == k * k * tt

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

    precision <- if (is.null(precision_mean)) {
      diag(1, k)
    } else if (precision_path) {
      matrix(precision_mean[(t - 1) * k * k + seq_len(k * k)], k)
    } else {
      matrix(precision_mean, k)
    }
    metric_u <- .symmetric_roots(precision)

    decomposition <- svd(metric_u[["root"]] %*% pi_mean %*% metric_w[["root"]], nu = rank, nv = rank)
    alpha_t <- metric_u[["inv_root"]] %*% decomposition[["u"]] %*%
      diag(decomposition[["d"]][seq_len(rank)], rank)
    beta_t <- metric_w[["inv_root"]] %*% decomposition[["v"]]

    a_point[pos_a + seq_len(n_alpha)] <- as.numeric(alpha_t)
    beta_point[pos_beta + seq_len(n_beta)] <- as.numeric(beta_t)
  }

  point[["a"]][["coeffs"]] <- coda::mcmc(matrix(a_point, nrow = 1), start = 1, end = 1, thin = 1)
  point[["beta"]][["coeffs"]] <- coda::mcmc(matrix(beta_point, nrow = 1), start = 1, end = 1, thin = 1)

  point
}

# The symmetric square root of a positive semi-definite matrix and the
# Moore-Penrose inverse of that root. Directions without variation -- a series
# in W that is a combination of others -- get a root of zero, and the fit in
# them is indifferent to Pi, so the pseudo-inverse leaves them out.
.symmetric_roots <- function(m) {
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  values <- pmax(e[["values"]], 0)
  keep <- values > max(values) * 1e-12
  vectors <- e[["vectors"]]
  list("root" = vectors %*% diag(sqrt(values), length(values)) %*% t(vectors),
       "inv_root" = vectors[, keep, drop = FALSE] %*%
         diag(1 / sqrt(values[keep]), sum(keep)) %*% t(vectors[, keep, drop = FALSE]))
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
