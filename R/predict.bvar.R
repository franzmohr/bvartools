#' Predict Method for Objects of Class bvar
#' 
#' Forecasting a Bayesian VAR object of class "bvar" with credible bands.
#' 
#' @param object an object of class "bvar", usually, a result of a call to `bvar` or `bvec_to_bvar`.
#' @param n.ahead number of steps ahead at which to predict.
#' @param new_x a matrix of new non-deterministic, exogenous variables. Must have `n.ahead` rows.
#' @param new_D a matrix of new deterministic variables. Must have `n.ahead` rows.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#' 
#' @details For the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + \sum_{i = 0}^{s} B_{i} x_{t-i} + C D_t + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)} the function produces `n.ahead` forecasts.
#' 
#' @return A time-series object of class "irf.bvar".
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' @export
#' @rdname bvar
predict.bvar <- function(object, ..., n.ahead = 10, new_y = NULL, new_x = NULL, new_D = NULL, ci = .95) {
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
      #stop("For the inclusion of exogenous please specify the 'new_x' argument.")
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
      #stop("For the inclusion of deterministic terms please specify the 'new_D' argument.")
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
  if (is.null(new_y)) {
    pred[pos_y, 1] <- object$y[, t:(t - p + 1)]
  } else {
    pred[pos_y, 1] <- new_y
  }
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
      if (struct) {
        A0_i <- solve(matrix(object$A0[draw, ], k))
      }
      temp <- svd(matrix(object$Sigma[draw, ], k))
      u <- temp$u %*% diag(sqrt(temp$d), k) %*% t(temp$v) %*% stats::rnorm(k)
      pred[1:k, i + 1] <- A0_i %*% matrix(A[draw, ], k) %*% pred[, i] + u
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
  
  result <- list("y" = object$y,
                 "fcst" = result)
  
  class(result) <- append("bvarprd", class(result))
  return(result)
}