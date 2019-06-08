#' Predict Method for Objects of Class bvar
#' 
#' Forecasting a Bayesian VAR object of class \code{"bvar"} with credible bands.
#' 
#' @param object an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param n.ahead number of steps ahead at which to predict.
#' @param new_x a matrix of new non-deterministic, exogenous variables. Must have \code{n.ahead} rows.
#' @param new_D a matrix of new deterministic variables. Must have \code{n.ahead} rows.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#' 
#' @details For the VAR model
#' \deqn{y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + \sum_{i = 0}^{s} B_{i} x_{t-i} + C D_t + A_0^{-1} u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)} the function produces \code{n.ahead} forecasts.
#' 
#' @return A time-series object of class \code{"bvarprd"}.
#' 
#' @examples
#' data("e1")
#' e1 <- diff(log(e1))
#' 
#' data <- gen_var(e1, p = 2, deterministic = "const")
#' 
#' y <- data$Y[, 1:73]
#' x <- data$Z[, 1:73]
#' 
#' set.seed(1234567)
#' 
#' iter <- 3000 # Number of iterations of the Gibbs sampler
#' burnin <- 1000 # Number of burn-in draws
#' store <- iter - burnin
#' 
#' t <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' m <- k * nrow(x) # Number of estimated coefficients
#' 
#' # Set (uninformative) priors
#' a_mu_prior <- matrix(0, m) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, m) # Inverse of the prior covariance matrix
#' 
#' u_sigma_df_prior <- 0 # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
#' u_sigma_df_post <- t + u_sigma_df_prior # Posterior degrees of freedom
#' 
#' # Initial values
#' u_sigma_i <- diag(.00001, k)
#' u_sigma <- solve(u_sigma_i)
#' 
#' # Data containers for posterior draws
#' draws_a <- matrix(NA, m, store)
#' draws_sigma <- matrix(NA, k^2, store)
#' 
#' # Start Gibbs sampler
#' for (draw in 1:iter) {
#'   # Draw conditional mean parameters
#'   a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
#' 
#' # Draw variance-covariance matrix
#' u <- y - matrix(a, k) %*% x # Obtain residuals
#' u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
#' u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#' u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
#' 
#' # Store draws
#' if (draw > burnin) {
#'   draws_a[, draw - burnin] <- a
#'   draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#' 
#' # Generate bvar object
#' bvar_est <- bvar(y = y, x = x, A = draws_a[1:18,],
#'                  C = draws_a[19:21, ], Sigma = draws_sigma)
#' 
#' # Generate forecasts
#' bvar_pred <- predict(bvar_est, n.ahead = 10, new_D = rep(1, 10))
#' 
#' # Plot forecasts
#' plot(bvar_pred)
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
#' @rdname bvar
predict.bvar <- function(object, ..., n.ahead = 10, new_x = NULL, new_D = NULL, ci = .95) {
  k <- object$specifications$dims["K"]
  t <- ncol(object$y)
  A <- object$A
  p <- object$specifications$lags["p"]
  tot <- k * p
  
  if (!is.null(object$B)) {
    A <- cbind(A, object$B)
    m <- ncol(object$B) / k
    tot <- tot + m
    if (is.null(new_x)) {
      new_x <- matrix(0, n.ahead, m)
    }
    if (NROW(new_x) != n.ahead) {
      stop("Length of argument 'new_x' must be equal to 'n.ahead'.")
    }
  }
  if (!is.null(object$C)) {
    A <- cbind(A, object$C)
    n <- ncol(object$C) / k
    tot <- tot + n
    if (is.null(new_D)) {
      new_D <- matrix(0, n.ahead, n)
    }
    if (NROW(new_D) != n.ahead) {
      stop("Length of argument 'new_D' must be equal to 'n.ahead'.")
    }
  }
  
  if (is.null(object$A0)) {
    struct <- FALSE
  } else {
    struct <- TRUE
  }
  
  pos_y <- 1:(k * p)
  pred <- matrix(NA, tot, n.ahead + 1)
  pred[pos_y, 1] <- object$y[, t:(t - p + 1)]
  if (!is.null(new_x)) {
    pos_x <- k * p + 1:m
    pred[pos_x, ] <- cbind(object$x[pos_x, t], t(new_x))
  }
  if (!is.null(object$C)) {
    pos_d <- (tot - n + 1):tot
    pred[pos_d, ] <- cbind(object$x[pos_d, t], t(new_D))
  }
  
  A0_i <- diag(1, k)
  draws <- nrow(A)
  result <- array(NA, dim = c(k, n.ahead, draws))
  
  for (draw in 1:draws) {
    for (i in 1:n.ahead) {
      temp <- eigen(matrix(object$Sigma[draw, ], k))
      u <- temp$vectors %*% diag(sqrt(temp$values), k) %*% t(temp$vectors) %*% stats::rnorm(k)
      if (struct) {
        A0_i <- solve(matrix(object$A0[draw, ], k))
      }
      pred[1:k, i + 1] <- matrix(A[draw, ], k) %*% pred[, i] + A0_i %*% u
      for (j in 1:(p - 1)) {
        pred[j * k + 1:k, i + 1] <- pred[(j - 1) * k + 1:k, i]
      }
    }
    result[,, draw] <- pred[1:k, -1]
  }
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(result, c(2, 1) , stats::quantile, probs = c(ci_low, .5, ci_high)) 
  result <- c()
  for (i in 1:k) {
    result <- c(result, list(stats::ts(t(temp[,, i]))))
  }
  names(result) <- dimnames(object$y)[[1]]
  
  if (!is.null(attr(object$y, "ts_info"))) {
    ts_info <- attr(object$y, "ts_info")
    object$y <- stats::ts(t(object$y), start = ts_info[1], frequency = ts_info[3])
    attr(object$y, "ts_info") <- NULL
    
    ts_temp <- stats::ts(0:n.ahead, start = ts_info[2], frequency = ts_info[3])
    ts_temp <- stats::time(ts_temp)[-1]
    for (i in 1:k) {
      stats::tsp(result[[i]]) <- c(ts_temp[1], ts_temp[length(ts_temp)], ts_info[3])
    }
  } else {
    object$y <- stats::ts(t(object$y))
    ts_temp <- stats::tsp(object$y)
    for (i in 1:k) {
      result[[i]] <- stats::ts(result[[i]], start = ts_temp[2] + 1, frequency = 1)
    }
  }
  
  result <- list("y" = object$y,
                 "fcst" = result)
  
  class(result) <- append("bvarprd", class(result))
  return(result)
}