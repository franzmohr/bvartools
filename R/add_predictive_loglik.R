#' Add Predictive Log-Likelihood
#'
#' Adds the draws of the one-step-ahead log predictive density to the windows of
#' an expanding window exercise, which \code{\link{selection_criteria}} sums to
#' the log predictive likelihood.
#'
#' @param object an object of class \code{"expandingwindow"} whose windows hold
#' posterior draws, or a \code{"modellist"} of such objects, as returned by
#' \code{\link{use_expanding_window}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The windows of an expanding window exercise grow by one period at a time, so
#' window \eqn{i + 1} holds exactly one observation that window \eqn{i} has not
#' seen. For every window but the last, the function evaluates the density of
#' that observation, \eqn{\Delta y_{t}}, given the regressors of period \eqn{t},
#' at each posterior draw of the window, and stores the draws in element
#' \code{predictive} of the window, together with the period \eqn{t}. The log of
#' their mean is the log predictive density \eqn{\ln p(\Delta y_t | y_{t-1},
#' \ldots, y_1)}, and its sum over the windows is the log predictive likelihood
#' of Geweke and Amisano (2011), which Koop, León-González and Strachan (2011)
#' use to choose between time varying cointegration models and their ranks.
#'
#' A draw of a model with constant coefficients and a constant error covariance
#' describes period \eqn{t} as it describes the sample. Everything that follows
#' a state equation is carried one period forward with that equation first:
#' \itemize{
#'  \item time varying coefficients by one step of their random walk with the
#' drawn state variance, and a coefficient that variable selection excluded
#' stays at zero;
#'  \item a time varying cointegration space by
#' \eqn{\beta_t = \rho (I_r \otimes P_\tau) \beta_{t-1} + \eta_t},
#' \eqn{\eta_t \sim N(0, I)}, with the drawn \eqn{\rho} or the one of the prior,
#' where \eqn{P_\tau} is the identity unless the prior centres the space;
#'  \item time varying error covariances by one step of their random walk;
#'  \item stochastic volatilities by one step of the random walk of the log
#' variances. Its state variance is taken from \code{posterior$u_sigma_inv$sigma}
#' if the sampler stored it, and is otherwise drawn from its conditional
#' posterior given the drawn path of the log variances and the prior in
#' \code{priors$u_sigma}.
#' }
#'
#' Unlike the pointwise log-likelihood of \code{\link{add_posterior_loglik}},
#' whose states have seen the observation they are evaluated at, the predictive
#' density conditions only on the data before period \eqn{t}. This is what makes
#' it the criterion for models whose coefficients or variances follow a state
#' equation, where leave-one-out importance sampling fails for exactly the
#' periods a path bends towards.
#'
#' The function is available for VEC models, which includes the rank zero
#' models in differences that ranks are compared with, and for VAR models, whose
#' regressors are the lags of the endogenous variables, the exogenous variables
#' and the deterministic terms of the period that is predicted. A VAR model is
#' the case of a VEC model without an error correction term, so the density is
#' the same expression with the cointegration block left out. The windows of a
#' VEC model must have been simulated on error correction terms that are neither
#' scaled nor centred, or put back with
#' \code{\link{rescale_error_correction}} first, and structural models are not
#' supported.
#'
#' The expression is the normal density of the observation, so the function is
#' available for the algorithms whose observation is normal given the parameters
#' of a draw, and refuses the others. The asymmetric Laplace algorithms of
#' quantile estimation, \code{"VarNormalAld"} and \code{"VarTvpAld"}, are
#' refused: the precision their samplers store is the one of the normal that
#' their scale mixture conditions on period by period, not the density of an
#' observation, whose mixing variable would have to be integrated out.
#'
#' The methods for a single fitted model, \code{\link{add_predictive_loglik.bvarmodel}}
#' and \code{\link{add_predictive_loglik.bvecmodel}}, score a forecast against
#' the observations its horizon realised instead. That is the same statistic over
#' a different set of periods -- one step ahead densities, each conditioning on
#' the realised history before it -- and it is BayesTS that computes it rather
#' than the R code here. It is stored in \code{posterior$forecast$loglik}, beside
#' the forecasts it scores, rather than in \code{predictive}.
#'
#' @return The object in \code{object}, with element \code{predictive} added to
#' every window but the last. It is a list with \code{loglik}, the draws of the
#' log predictive density, and \code{period}, the time of the predicted
#' observation.
#'
#' @seealso Methods for a single fitted model:
#' \code{\link{add_predictive_loglik.bvarmodel}},
#' \code{\link{add_predictive_loglik.bvecmodel}}.
#'
#' @references
#'
#' Geweke, J., & Amisano, G. (2011). Hierarchical Markov normal mixture models with
#' applications to financial asset returns. \emph{Journal of Applied Econometrics, 26}(1),
#' 1--29. \doi{10.1002/jae.1119}
#'
#' Koop, G., León-González, R., & Strachan, R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#'
#' @examples
#'
#' # Load data
#' data("e6")
#' e6 <- e6 * 100
#'
#' # Create model
#' model <- create_bvecmodel(e6, p = 2, r = 1, const = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 0.0001))
#'
#' # Estimate the model on the last four expanding windows
#' model <- use_expanding_window(model, start = c(1998, 1))
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' # One-step-ahead log predictive densities
#' model <- add_predictive_loglik(model)
#'
#' # The log predictive likelihood is criterion "LPL"
#' selection_criteria(model)
#'
#' @family model comparison
#' @export
add_predictive_loglik <- function(object, ...) {
  UseMethod("add_predictive_loglik")
}

#' @rdname add_predictive_loglik
#' @export
add_predictive_loglik.expandingwindow <- function(object, ...) {

  n_windows <- length(object)
  if (n_windows < 2) {
    stop("An expanding window needs at least two windows: the last one only ",
         "provides the observation the one before it predicts.")
  }

  for (i in seq_len(n_windows)) {
    if (!inherits(object[[i]], "bvecmodel") && !inherits(object[[i]], "bvarmodel")) {
      stop("Predictive log-likelihoods are available for VAR and VEC models only.")
    }
    .check_predictive_algorithm(object[[i]], i)
  }

  for (i in seq_len(n_windows - 1)) {
    current <- object[[i]]
    following <- object[[i + 1]]
    .check_predictive_window(current, following, i)

    newdata <- .predictive_observation(following, .predictive_rank(current))
    object[[i]][["predictive"]] <- list(
      loglik = .predictive_log_density(current, newdata),
      period = newdata[["period"]])
  }

  return(object)
}

#' @rdname add_predictive_loglik
#' @export
add_predictive_loglik.modellist <- function(object, ...) {

  orig_class <- class(object)
  object <- lapply(object, add_predictive_loglik, ...)
  class(object) <- orig_class

  return(object)
}

# The algorithms the density below is the density of: the ones whose observation
# is normal given the parameters of a draw. Everything else is refused rather
# than handed to an expression that does not describe it, which a whitelist does
# and a list of exceptions would not -- an algorithm added later is refused
# until its density is written down.
#
# The asymmetric Laplace models are the ones this keeps out today. They do store
# 'u_sigma_inv', the precision of the normal that the sampler's scale mixture
# conditions on period by period, so the expression below would run on them and
# return a number; but that number is a normal density at the mixing variables
# of the last period of the sample, not the predictive density of an asymmetric
# Laplace observation, whose mixing variable of the predicted period has to be
# integrated out.
# It governs both statistics: the expanding window density computed in R below
# and the score BayesTS takes of a forecast in add_predictive_loglik.bvarmodel().
# Upstream implements predictive_log_density() for exactly these thirteen
# algorithms and for no others, so the two lists are the same list.
.predictive_algorithms <- c("VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                            "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart",
                            "VecKlgs2010", "VecNormalGamma", "VecNormalStochvol",
                            "VecNormalWishart", "VecTvpGamma", "VecTvpStochvol",
                            "VecTvpWishart", "VarTvpDiscount", "VecTvpDiscount")

# Refuses a model whose algorithm the density is not the density of. 'i' names
# the window it is, where the caller is working through an expanding window
# exercise, and is left out where the model stands on its own.
.check_predictive_algorithm <- function(model, i = NULL) {

  algorithm <- model[["model"]][["algorithm"]]
  if (is.null(algorithm) || !algorithm %in% .predictive_algorithms) {
    stop("Predictive log-likelihoods are not available for the algorithm ",
         if (is.null(algorithm)) "" else paste0("'", algorithm, "' "),
         if (is.null(i)) "of this model" else paste0("of window ", i),
         ". They are available for the algorithms with a normal likelihood: ",
         paste(.predictive_algorithms, collapse = ", "), ".")
  }

  invisible(TRUE)
}

# The rank of a window. A VAR model has none, and enters the density as the case
# of rank zero.
.predictive_rank <- function(model) {
  rank <- model[["model"]][["rank"]]
  if (is.null(rank)) 0L else as.integer(rank)
}

# Refuses a pair of windows the density cannot be evaluated for, and says which.
.check_predictive_window <- function(current, following, i) {

  if (isTRUE(current[["model"]][["structural"]])) {
    stop("Predictive log-likelihoods are not available for structural models.")
  }
  if (is.null(current[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Window ", i, " holds no posterior draws. Use 'add_posterior_coefficients' first.")
  }
  transformed <- function(model) {
    w <- model[["data"]][["train"]][["w"]]
    !is.null(attr(w, "scale")) || !is.null(attr(w, "centre"))
  }
  # Only a VEC model has an error correction term to transform.
  if (inherits(current, "bvecmodel") && (transformed(current) || transformed(following))) {
    stop("The error correction terms of window ", i, " are scaled or centred. Use ",
         "'rescale_error_correction' before 'add_predictive_loglik'.")
  }

  y_current <- current[["data"]][["train"]][["y"]]
  y_following <- following[["data"]][["train"]][["y"]]
  tsp_current <- stats::tsp(y_current)
  tsp_following <- stats::tsp(y_following)
  if (NROW(y_following) != NROW(y_current) + 1 ||
      tsp_following[1] != tsp_current[1] ||
      abs(tsp_following[2] - tsp_current[2] - 1 / tsp_current[3]) > 1e-8) {
    stop("Window ", i + 1, " does not extend window ", i, " by exactly one ",
         "period. Use 'use_expanding_window' to create the windows.")
  }

  invisible(TRUE)
}

# The observation the following window adds: its last row of the response, of
# the error correction term and of the other regressors.
.predictive_observation <- function(model, rank) {

  y <- model[["data"]][["train"]][["y"]]
  tt <- NROW(y)
  x <- model[["data"]][["train"]][["x"]]
  w <- model[["data"]][["train"]][["w"]]

  list(y = as.numeric(y[tt, ]),
       w = if (rank > 0) as.numeric(w[tt, ]) else NULL,
       x = if (!is.null(x) && NCOL(x) > 0) as.numeric(x[tt, ]) else numeric(0),
       period = stats::time(y)[tt])
}

# Draws of the log density of one observation under the posterior of a VAR or a
# VEC model.
#
# 'newdata' holds the observation, y, and its regressors: w, the error
# correction term (NULL at rank zero), and x, the other regressors, both in the
# layout of the model's data$train. With 'innovations = FALSE' no state is
# carried forward, which makes the result the pointwise log-likelihood of the
# last period of the sample if 'newdata' is that period -- the identity the
# tests check the layouts against.
.predictive_log_density <- function(model, newdata, innovations = TRUE) {

  post <- model[["posterior"]]
  k <- NCOL(model[["data"]][["train"]][["y"]])
  kk <- k * k
  r <- .predictive_rank(model)
  tvp <- isTRUE(model[["model"]][["tvp"]])
  sv <- model[["model"]][["error"]] %in% c("sv", "sv+covar")
  tt <- NROW(model[["data"]][["train"]][["y"]])

  x <- newdata[["x"]]
  w <- newdata[["w"]]
  n_x <- length(x)
  k_w <- length(w)
  n_a <- k * (r + n_x)
  nd <- NROW(post[["u_sigma_inv"]][["coeffs"]])

  # The draws of the last period of a path, or the draws themselves where the
  # block does not move over time.
  last <- function(draws, size) {
    draws <- as.matrix(draws)
    if (ncol(draws) == size) {
      return(draws)
    }
    draws[, ncol(draws) - size + seq_len(size), drop = FALSE]
  }
  step <- function(sd) {
    matrix(stats::rnorm(length(sd)), NROW(sd)) * sd
  }

  # Coefficients
  if (n_a > 0) {
    a <- last(post[["a"]][["coeffs"]], n_a)
    if (tvp && innovations && !is.null(post[["a"]][["sigma"]])) {
      a_step <- step(sqrt(as.matrix(post[["a"]][["sigma"]])))
      # An excluded coefficient is zero in every period and stays there.
      if (!is.null(post[["a"]][["lambda"]])) {
        a_step <- a_step * as.matrix(post[["a"]][["lambda"]])
      }
      a <- a + a_step
    }
  }

  # Cointegration space
  if (r > 0) {
    n_beta <- k_w * r
    beta <- last(post[["beta"]][["coeffs"]], n_beta)
    if (tvp && innovations) {
      rho <- if (!is.null(post[["beta"]][["rho"]])) {
        as.numeric(post[["beta"]][["rho"]])
      } else {
        rep(model[["priors"]][["beta"]][["rho"]], nd)
      }
      p_tau <- model[["priors"]][["beta"]][["p_tau"]]
      if (!is.null(p_tau)) {
        beta <- beta %*% t(kronecker(diag(1, r), p_tau))
      }
      beta <- rho * beta + matrix(stats::rnorm(nd * n_beta), nd)
    }
  }

  # Error precision
  if (!is.null(post[["u_omega_inv"]][["coeffs"]])) {
    omega_inv <- last(post[["u_omega_inv"]][["coeffs"]], k)
    if (sv && innovations) {
      h_variance <- .sv_state_variance(model, nd, k, tt)
      omega_inv <- omega_inv * exp(-step(sqrt(h_variance)))
    }

    psi <- if (!is.null(post[["psi"]][["coeffs"]])) {
      last(post[["psi"]][["coeffs"]], kk)
    } else {
      matrix(as.numeric(diag(1, k)), nd, kk, byrow = TRUE)
    }
    # The free elements of Psi, below the diagonal, row by row, which is the
    # order of the draws of their state variances.
    if (tvp && innovations && k > 1 && !is.null(post[["psi"]][["sigma"]])) {
      free <- unlist(lapply(2:k, function(i) (seq_len(i - 1) - 1) * k + i))
      psi_step <- step(sqrt(as.matrix(post[["psi"]][["sigma"]])))
      if (!is.null(post[["psi"]][["lambda"]])) {
        psi_step <- psi_step * as.matrix(post[["psi"]][["lambda"]])[, free, drop = FALSE]
      }
      psi[, free] <- psi[, free] + psi_step
    }

    precision <- function(d) {
      p <- matrix(psi[d, ], k)
      crossprod(p, omega_inv[d, ] * p)
    }
  } else {
    u_sigma_inv <- last(post[["u_sigma_inv"]][["coeffs"]], kk)
    precision <- function(d) {
      matrix(u_sigma_inv[d, ], k)
    }
  }

  result <- vapply(seq_len(nd), function(d) {
    mu <- numeric(k)
    if (n_a > 0) {
      coefficients <- matrix(a[d, ], k)
      if (r > 0) {
        mu <- mu + coefficients[, seq_len(r), drop = FALSE] %*%
          crossprod(matrix(beta[d, ], k_w), w)
      }
      if (n_x > 0) {
        mu <- mu + coefficients[, r + seq_len(n_x), drop = FALSE] %*% x
      }
    }
    u <- newdata[["y"]] - as.numeric(mu)
    p <- precision(d)
    -k / 2 * log(2 * pi) +
      as.numeric(determinant(p, logarithm = TRUE)[["modulus"]]) / 2 -
      as.numeric(crossprod(u, p %*% u)) / 2
  }, numeric(1))

  return(result)
}

# Draws of the variance of the random walk of the log variances, draws x k.
#
# Taken from the posterior if the sampler stored them. Otherwise each draw is
# made from its conditional posterior given the drawn path of the log variances,
# an inverse gamma with the prior of priors$u_sigma, which conditions on the
# path from its second period on and so leaves out the state before the sample.
.sv_state_variance <- function(model, nd, k, tt) {

  stored <- model[["posterior"]][["u_sigma_inv"]][["sigma"]]
  if (!is.null(stored)) {
    return(as.matrix(stored))
  }

  if (tt < 2) {
    stop("The variance of the log volatilities cannot be drawn from a sample of one period.")
  }
  h <- array(-log(as.matrix(model[["posterior"]][["u_omega_inv"]][["coeffs"]])),
             c(nd, k, tt))
  squares <- apply((h[, , -1, drop = FALSE] - h[, , -tt, drop = FALSE])^2, c(1, 2), sum)
  shape <- rep(as.numeric(model[["priors"]][["u_sigma"]][["shape"]]), length.out = k)
  rate <- rep(as.numeric(model[["priors"]][["u_sigma"]][["rate"]]), length.out = k)
  precision <- stats::rgamma(nd * k,
                             shape = rep(shape + (tt - 1) / 2, each = nd),
                             rate = rep(rate, each = nd) + as.numeric(squares) / 2)

  return(matrix(1 / precision, nd, k))
}

# Log of the mean of exp(x), and its numerical standard error by the delta
# method, with batch means for the autocorrelation of the chain.
.log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# The 'LPL' entry of a set of criteria, built from a sequence of one step ahead
# predictive densities.
#
# 'densities' is a list of draws of the log predictive density, one element per
# period, and 'periods' the times of the observations they score. The log of the
# mean of exp() of each is that period's log predictive density, and their sum
# is the log predictive likelihood.
#
# Both sources reach this. An expanding window exercise has one density per
# window, taken by add_predictive_loglik() from the draws of that window; a
# single model has one per horizon of its forecast, taken by BayesTS against the
# values in data$test$y and read back from posterior$forecast$loglik. They are
# the same quantity computed two ways, and going through one function is what
# makes them agree by construction rather than by coincidence.
#
# The band treats the terms as independent, which they are not -- consecutive
# windows share all but one observation, and consecutive horizons share a
# forecast origin -- so it is a rough guide to how much of the total is sampling
# noise rather than an interval to test with. The numerical standard error
# beside it is a different thing: what the length of the chain contributes, which
# more draws would shrink and more data would not.
.lpl_entry <- function(densities, periods, ci_low, ci_high) {

  terms <- data.frame(
    period = as.numeric(periods),
    lpd = vapply(densities, .log_mean_exp, numeric(1)),
    nse = vapply(densities, .nse_log_mean_exp, numeric(1)))

  n_terms <- nrow(terms)
  lpl <- sum(terms[["lpd"]])
  se <- if (n_terms > 1) sqrt(n_terms * stats::var(terms[["lpd"]])) else NA_real_
  z <- stats::qnorm(ci_high)

  entry <- data.frame(mean = lpl, median = NA_real_,
                      qlower = lpl - z * se, qupper = lpl + z * se)
  attr(entry, "terms") <- terms
  attr(entry, "nse") <- sqrt(sum(terms[["nse"]]^2))

  return(entry)
}


# The draws of the log predictive density a model carries, one element per
# scored period, with the times of the periods they score. NULL where the model
# has not been scored.
#
# posterior$forecast$loglik is draws by scored periods, written by BayesTS
# against data$test$y. The periods are counted on from the end of the estimation
# sample, and fall back on the horizon number where the sample carries no time.
.model_predictive_densities <- function(object) {

  score <- object[["posterior"]][["forecast"]][["loglik"]]
  if (is.null(score)) {
    return(NULL)
  }

  score <- as.matrix(score)
  n <- ncol(score)
  tsp_train <- stats::tsp(object[["data"]][["train"]][["y"]])
  periods <- if (is.null(tsp_train)) {
    seq_len(n)
  } else {
    tsp_train[2] + seq_len(n) / tsp_train[3]
  }

  return(list(densities = lapply(seq_len(n), function(i) score[, i]),
              periods = periods))
}


.nse_log_mean_exp <- function(x, batches = 20) {
  batches <- min(batches, floor(length(x) / 2))
  if (batches < 2) {
    return(NA_real_)
  }
  e <- exp(x - max(x))
  means <- vapply(split(e, cut(seq_along(e), batches, labels = FALSE)), mean, numeric(1))
  sqrt(stats::var(means) / batches) / mean(e)
}
