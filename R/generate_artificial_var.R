#' Generate Artificial VAR Data
#' 
#' Generates an artificial data set for algorithm testing.
#' 
#' @param nobs number of generated observations. Defaults to 100.
#' @param k number of endogenous variables. Defaults to 3.
#' @param p number of lags of the VAR model. Defaults to 2.
#' @param const logical specifying if the model should include an intercept term.
#' Defaults to \code{FALSE}.
#' @param a_range numeric vector with two elements containing the minimum and
#' maximum value of the coefficients. Defaults to \code{c(-0.5, 0.5)}.
#' @param a_zeros numeric between 0 and 1 indicating the share of parameters,
#' which should be set to zero. Default is \code{0.5}.
#' maximum value of the coefficients. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_const numeric vector with two elements containing the minimum and
#' maximum value of the coefficients of intercept terms. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_variance numeric vector with two elements containing the minimum and
#' maximum value of the error variances. Defaults to \code{c(1, 1)}.
#' @param range_psi numeric vector with two elements containing the minimum and
#' maximum value of the error covariances. Defaults to \code{c(0, 0)}.
# @param range_variance_state numeric vector with two elements containing the minimum and
# maximum value of the error variances of the state equation. Only used for
# state-space models.
#' @param presample numeric specifying the number of presample values used to
#' initialize the time series. Defaults to zero.
#' 
#' @details The function produces artificial observations for a vector
#' autoregressive (VAR) model:
#' \deqn{y_t = c + \sum_{i=1}^{p} A_i y_{t - i} + u_t,}
#' where
#' \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{c} is a \eqn{K \times 1} vector of intercept terms.
#' \eqn{p} is the lag order of endogenous variables.
#' 
#' As in Primiceri (2005) \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Psi \Omega_u \Psi^{\prime})},
#' where \eqn{\Psi} is a lower triangular matrix with ones on the main diagonal
#' and \eqn{\Omega} is a diagonal matrix with error variances.
#' 
#' @references
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Primiceri, G. E. (2005). Time varying structural vector autoregressions and
#' monetary policy. \emph{The Review of Economic Studies, 72}(3), 821--852.
#' \doi{10.1111/j.1467-937X.2005.00353.x}
#' 
#' @returns A list containing the artificial time series and a list with the
#' true parameters.
#' 
#' @examples
#' 
#' # Set seed of RNG
#' set.seed(1)
#' 
#' # Time series without intercept terms
#' dt <- generate_artificial_var(nobs = 200, k = 3)
#' 
#' # Time series with all intercept terms equal to 5
#' dt <- generate_artificial_var(nobs = 200, k = 3
#'                               const = TRUE, range_const = c(5, 5))
#' 
#' 
#' @export
generate_artificial_var <- function(nobs = 100, k = 3, p = 2,
                                    const = FALSE,
                                    a_range = c(-0.5, 0.5),
                                    a_zeros = .5,
                                    range_const = NULL,
                                    range_variance = c(1, 1),
                                    range_psi = c(0, 0),
                                    presample = 0) {
  
  # Basic checks
  a_range <- c(min(a_range), max(a_range))
  if (a_zeros < 0 | a_zeros > 1) {
    stop("Argument 'a_zeros' must be between 0 and 1.")
  }
  a_prob_0 <- a_zeros
  a_prob_1 <- 1 - a_zeros
  if (const) {
    if (is.null(range_const)) {
      stop("If argument 'const' is TRUE, argument 'range_const' must be specified.")
    }
    range_const <- c(min(range_const), max(range_const)) 
  }
  range_variance <- c(min(range_variance), max(range_variance))
  range_psi <- c(min(range_psi), max(range_psi))
  #if (!is.null(range_variance_state)) {
  #  range_variance_state <- c(min(range_variance_state), max(range_variance_state)) 
  #}
  
  if (p > 0 | const) {
    # Create a
    n_a <- k * k * p
    
    not_stable <- TRUE
    while (not_stable) {
      a <- NULL
      # Generate candidate
      if (p > 0) {
        
        a_new <- stats::runif(n_a, a_range[1], a_range[2])
        
        a_pos_zeros <- which(sample(c(0, 1), n_a, prob = c(a_prob_0, a_prob_1), replace = TRUE) == 0)
        if (length(a_pos_zeros) > 0) {
          a_new[a_pos_zeros] <- 0 
        }
        
        a <- cbind(a, round(matrix(a_new, k), 2)) # Generate candidate
      }
      # Check stationarity
      if (nrow(a) < ncol(a)) {
        a_plus <- matrix(0, ncol(a) - nrow(a), ncol(a))
        a_plus[1:(ncol(a) - nrow(a)), 1:(ncol(a) - nrow(a))] <- diag(1, ncol(a) - nrow(a))
        not_stable <- any(abs(eigen(rbind(a, a_plus))$values) >= 1)
      } else {
        not_stable <- any(abs(eigen(a)$values) >= 1) 
      }
      if (const) {
        a <- cbind(a, round(matrix(stats::runif(k, range_const[1], range_const[2]), k), 2))
        pos_const <- ncol(a) # Save position of intercept coefficients for later
      }
    }
  } else {
    a <- NULL
  }
  
  # Create sigma
  omega <- diag(round(stats::runif(k, range_variance[1], range_variance[2]), 2), k)
  
  if (k > 1) {
    psi <- diag(1, k)
    psi[lower.tri(psi)] <- round(stats::runif(k * (k - 1) / 2, range_psi[1], range_psi[2]), 2)
    psi_inv <- solve(psi)
    sigma <- psi_inv %*% omega %*% t(psi_inv)
  } else {
    psi <- NULL
    sigma <- omega
  }
  
  # Generate time series
  y <- matrix(0, k, nobs + p + presample)
  if (p > 0) {
    if (const) {
      init_cond <- diag(1, k)
      for (i in 1:p) {
        init_cond <- init_cond - a[, k * (i - 1) + 1:k]
      }
      y[, 1:p] <- solve(init_cond) %*% a[, pos_const]
    }
    y[, 1:p] <- y[, 1:p] + stats::rnorm(k * p)
  }
  
  for (i in p + 1:(nobs + presample)) {
    if (const) {
      y[, i] <- a[, pos_const]
    }
    if (p > 0) {
      for (j in 1:p) {
        y[, i] <- y[, i] + a[, k * (j - 1) + 1:k] %*% y[, i - j] 
      } 
    }
    y[, i] <- y[, i] + chol(sigma) %*% matrix(stats::rnorm(k))
  }
  
  if (p > 0) {
    y <- y[, -(1:p)] # Drop first p observations 
  }
  if (presample > 0) {
    y <- y[, -(1:presample)] # Drop first p observations 
  }
  
  result <- list(data = stats::ts(t(y)),
                 params = list(a_coef = a,
                               psi_coef = psi,
                               u_omega = omega,
                               u_sigma = sigma))
  
  # Generate artificial names
  names_series <- paste0("var", 1:k)
  dimnames(result[["data"]]) <- list(NULL, names_series)
  
  return(result)
}