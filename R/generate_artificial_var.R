#' Generate Artificial VAR Data
#'
#' Generates an artificial data set from a vector autoregressive model with
#' constant or time varying parameters and constant or stochastic volatility,
#' for testing the algorithms of the package.
#'
#' @param nobs number of generated observations. Defaults to 100.
#' @param k number of endogenous variables. Defaults to 3.
#' @param p number of lags of the VAR model. Can be zero. Defaults to 2.
#' @param deterministic a character specifying which deterministic terms the
#' model contains: \code{"none"} (default), \code{"const"} for an intercept,
#' \code{"trend"} for a linear trend, and \code{"both"} for an intercept with a
#' linear trend. The terms are the same as in \code{\link{create_bvarmodel}}.
#' @param structural logical specifying whether a structural VAR model with a lower
#' triangular matrix \eqn{A_0} of contemporaneous coefficients is generated.
#' Defaults to \code{FALSE}. See 'Details'.
#' @param tvp logical specifying whether the coefficients, \eqn{A_0} and the
#' elements of \eqn{\Psi} follow random walks. Defaults to \code{FALSE}. See 'Details'.
#' @param sv logical specifying whether the error variances follow a stochastic
#' volatility process. Defaults to \code{FALSE}. See 'Details'.
#' @param range_a numeric vector with two elements containing the minimum and
#' maximum value of the coefficients of lagged endogenous variables. Defaults
#' to \code{c(-0.5, 0.5)}.
#' @param a_zeros numeric between 0 and 1 indicating the share of coefficients
#' of lagged endogenous variables, which should be set to zero. Default is \code{0.5}.
#' @param range_a0 numeric vector with two elements containing the minimum and
#' maximum value of the free elements of \eqn{A_0}. Only used if
#' \code{structural = TRUE}. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_const numeric vector with two elements containing the minimum and
#' maximum value of the intercept terms. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_trend numeric vector with two elements containing the minimum and
#' maximum value of the coefficients of the linear trend. Defaults to \code{c(-0.01, 0.01)}.
#' @param range_variance numeric vector with two elements containing the minimum and
#' maximum value of the error variances. For \code{sv = TRUE} these are the
#' variances in the first period. Defaults to \code{c(1, 1)}.
#' @param range_psi numeric vector with two elements containing the minimum and
#' maximum value of the free elements of \eqn{\Psi}. Defaults to \code{c(0, 0)},
#' which gives uncorrelated errors. Must be \code{c(0, 0)} for structural models.
#' @param range_variance_state numeric vector with two elements containing the minimum and
#' maximum value of the variances of the state equations of the coefficients, of
#' \eqn{A_0} and of \eqn{\Psi}. Only used if \code{tvp = TRUE}. Defaults to
#' \code{c(0.0001, 0.0001)}.
#' @param range_variance_sv numeric vector with two elements containing the minimum and
#' maximum value of the variances of the state equations of the log-volatilities.
#' Only used if \code{sv = TRUE}. Defaults to \code{c(0.01, 0.01)}.
#' @param stable logical specifying whether the coefficients of lagged endogenous
#' variables are restricted to a stable VAR process. Defaults to \code{TRUE}.
#' See 'Details'.
#' @param presample numeric specifying the number of observations, which are
#' generated before the first returned observation and then discarded, so that
#' the series do not depend on their initial values. Defaults to 100.
#' @param level numeric vector with one or \code{k} elements, which are added to
#' the generated series to shift their levels. Defaults to 0. A non-zero level
#' requires an intercept. See 'Details'.
#'
#' @details The function produces artificial observations for a vector
#' autoregressive (VAR) model:
#' \deqn{A_{0t} y_t = \sum_{i=1}^{p} A_{it} y_{t - i} + C_t d_t + u_t,}
#' where
#' \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_{0t}} is a \eqn{K \times K} matrix of contemporaneous coefficients,
#' \eqn{A_{it}} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{d_t} is a vector of deterministic terms and \eqn{C_t} its coefficient
#' matrix. \eqn{p} is the lag order of endogenous variables.
#'
#' As in Primiceri (2005) \eqn{u_t} is an error term with
#' \eqn{u_t \sim N(0, \Sigma_t)} and \eqn{\Psi_t \Sigma_t \Psi_t^{\prime} = \Omega_t},
#' where \eqn{\Psi_t} is a lower triangular matrix with ones on the main diagonal
#' and \eqn{\Omega_t} is a diagonal matrix with error variances
#' \eqn{\omega_{1t}, \dots, \omega_{Kt}}.
#'
#' Unless \code{structural = TRUE}, \eqn{A_{0t}} is an identity matrix. Otherwise it
#' is a lower triangular matrix with ones on its main diagonal, whose free elements
#' are drawn from \code{range_a0}, and the errors of the structural equations are
#' uncorrelated, so that \eqn{\Psi_t} is an identity matrix and
#' \eqn{\Sigma_t = \Omega_t}. This is the model that \code{\link{create_bvarmodel}}
#' produces with \code{structural = TRUE}. The coefficients \eqn{A_{it}} and \eqn{C_t}
#' are those of the structural equations, and the reduced form of the model has the
#' coefficients \eqn{A_{0t}^{-1} A_{it}} and the error covariance matrix
#' \eqn{A_{0t}^{-1} \Omega_t A_{0t}^{-1 \prime}}.
#'
#' The linear trend takes the value 1 in period \eqn{p + 1} of the returned series,
#' which is the first period that \code{\link{create_bvarmodel}} uses for estimation
#' with lag order \eqn{p}. The true coefficients of the trend are therefore those
#' that model estimates with the same lag order refer to.
#'
#' If \code{tvp = FALSE} and \code{sv = FALSE}, all parameters are constant over time.
#' If \code{tvp = TRUE}, the coefficients \eqn{a_t = vec(A_{1t}, \dots, A_{pt}, C_t)},
#' the free elements of \eqn{A_{0t}} and the free elements \eqn{\psi_t} of
#' \eqn{\Psi_t} follow random walks
#' \deqn{a_t = a_{t-1} + v_t, \quad v_t \sim N(0, Q_a), \qquad
#' \psi_t = \psi_{t-1} + w_t, \quad w_t \sim N(0, Q_\psi),}
#' with diagonal \eqn{Q_a} and \eqn{Q_\psi}, whose elements are drawn from
#' \code{range_variance_state}, and likewise for \eqn{A_{0t}}. Coefficients that are
#' zero in the first period, such as those set to zero by \code{a_zeros}, have a state
#' variance of zero and stay zero. With the default \code{range_psi = c(0, 0)} the
#' errors therefore remain uncorrelated.
#' If \code{sv = TRUE}, the log-volatilities follow random walks
#' \deqn{\ln \omega_{it} = \ln \omega_{i,t-1} + \eta_{it}, \quad \eta_{it} \sim N(0, \sigma^2_i),}
#' with \eqn{\sigma^2_i} drawn from \code{range_variance_sv}. These are the state
#' equations of the models that \code{\link{create_bvarmodel}} produces with
#' \code{tvp = TRUE} and \code{error = "sv"} or \code{"sv+covar"}.
#'
#' The parameters of the first period are drawn from the ranges above. During the
#' presample they are held at these values, so that the returned paths start at
#' them and the presample only removes the influence of the initial values of the
#' series.
#'
#' If \code{stable = TRUE}, the coefficients of lagged endogenous variables are
#' drawn until the companion matrix of the reduced form of the process has no
#' eigenvalue on or outside the unit circle. For \code{tvp = TRUE} this holds in
#' every period: an innovation of the coefficients, which would make the process
#' unstable, is drawn again, and if no stable innovation is found in 100 attempts,
#' the coefficients keep the values of the previous period. This restricts the
#' random walks to the stable region as in Cogley and Sargent (2005). Set
#' \code{stable = FALSE} to generate integrated or explosive series, usually together
#' with \code{presample = 0}. Cointegrated series are generated by
#' \code{\link{generate_artificial_vec}}.
#'
#' Argument \code{level} shifts the generated series by a vector \eqn{m}, which gives
#' them high levels without changing their dynamics. The model becomes
#' \eqn{A_{0t} (y_t - m) = \sum_{i=1}^{p} A_{it} (y_{t - i} - m) + C_t d_t + u_t},
#' so that the intercept absorbs the shift: the returned intercept is
#' \eqn{c_t + (A_{0t} - \sum_{i=1}^{p} A_{it}) m}, where \eqn{c_t} is drawn from
#' \code{range_const}. Without an intercept this term would not be part of the model,
#' so that a non-zero level requires \code{deterministic = "const"} or \code{"both"}.
#'
#' @references
#' Cogley, T., & Sargent, T. J. (2005). Drifts and volatilities: Monetary policies
#' and outcomes in the post WWII US. \emph{Review of Economic Dynamics, 8}(2), 262--302.
#' \doi{10.1016/j.red.2004.10.009}
#'
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' Primiceri, G. E. (2005). Time varying structural vector autoregressions and
#' monetary policy. \emph{The Review of Economic Studies, 72}(3), 821--852.
#' \doi{10.1111/j.1467-937X.2005.00353.x}
#'
#' @return A list with the elements
#' \describe{
#'   \item{\code{data}}{a \eqn{T \times K} time-series object of the artificial
#'   series, named \code{var1}, \code{var2} etc.}
#'   \item{\code{params}}{a list of the true parameters with the elements
#'   \describe{
#'     \item{\code{a_coef}}{the \eqn{K \times M} coefficient matrix
#'     \eqn{(A_1, \dots, A_p, C)}, with the columns in the order of the regressors of
#'     \code{\link{create_bvarmodel}}, or \code{NULL} if the model has no regressors.
#'     For \code{tvp = TRUE} a \eqn{K \times M \times T} array.}
#'     \item{\code{a0_coef}}{for \code{structural = TRUE}, the \eqn{K \times K} matrix
#'     \eqn{A_0}, or for \code{tvp = TRUE} a \eqn{K \times K \times T} array.}
#'     \item{\code{psi_coef}}{the \eqn{K \times K} matrix \eqn{\Psi}, or for
#'     \code{tvp = TRUE} a \eqn{K \times K \times T} array. \code{NULL} if \eqn{K = 1}
#'     or \code{structural = TRUE}.}
#'     \item{\code{u_omega}}{the \eqn{K \times K} diagonal matrix of error variances
#'     \eqn{\Omega}, or for \code{sv = TRUE} a \eqn{K \times K \times T} array.}
#'     \item{\code{u_sigma}}{the \eqn{K \times K} covariance matrix \eqn{\Sigma} of
#'     \eqn{u_t}, or for \code{tvp = TRUE} or \code{sv = TRUE} a
#'     \eqn{K \times K \times T} array.}
#'     \item{\code{a_state_variance}}{for \code{tvp = TRUE}, the \eqn{K \times M} matrix
#'     of the variances of the state equations of the coefficients.}
#'     \item{\code{a0_state_variance}}{for \code{tvp = TRUE} and \code{structural = TRUE},
#'     the \eqn{K \times K} matrix of the variances of the state equations of \eqn{A_0}.}
#'     \item{\code{psi_state_variance}}{for \code{tvp = TRUE}, \eqn{K > 1} and
#'     \code{structural = FALSE}, the \eqn{K \times K} matrix of the variances of the
#'     state equations of \eqn{\Psi}.}
#'     \item{\code{u_state_variance}}{for \code{sv = TRUE}, the \eqn{K} variances of
#'     the state equations of the log-volatilities.}
#'   }}
#' }
#' For time varying parameters, \code{as.vector()} of an array gives the parameters
#' period by period, which is the order of the columns of the posterior draws of TVP
#' models.
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
#' dt <- generate_artificial_var(nobs = 200, k = 3, deterministic = "const",
#'                               range_const = c(5, 5))
#'
#' # Time series with high levels, shifted by 100, 50 and 20
#' dt <- generate_artificial_var(nobs = 200, k = 3, deterministic = "const",
#'                               level = c(100, 50, 20))
#' dt[["params"]][["a_coef"]][, "const"]
#'
#' # Time varying parameters and stochastic volatility with correlated errors
#' dt <- generate_artificial_var(nobs = 200, k = 2, p = 1, deterministic = "const",
#'                               tvp = TRUE, sv = TRUE, range_psi = c(-0.5, 0.5))
#'
#' # Path of the error variance of the first variable
#' plot(dt[["params"]][["u_omega"]][1, 1, ], type = "l")
#'
#' # Structural model
#' dt <- generate_artificial_var(nobs = 200, k = 3, p = 1, structural = TRUE)
#' dt[["params"]][["a0_coef"]]
#'
#' @family artificial data
#' @export
generate_artificial_var <- function(nobs = 100, k = 3, p = 2,
                                    deterministic = "none",
                                    structural = FALSE,
                                    tvp = FALSE,
                                    sv = FALSE,
                                    range_a = c(-0.5, 0.5),
                                    a_zeros = 0.5,
                                    range_a0 = c(-0.5, 0.5),
                                    range_const = c(-0.5, 0.5),
                                    range_trend = c(-0.01, 0.01),
                                    range_variance = c(1, 1),
                                    range_psi = c(0, 0),
                                    range_variance_state = c(0.0001, 0.0001),
                                    range_variance_sv = c(0.01, 0.01),
                                    stable = TRUE,
                                    presample = 100,
                                    level = 0) {

  # Basic checks
  .artificial_count(nobs, "nobs", 1)
  .artificial_count(k, "k", 1)
  .artificial_count(p, "p", 0)
  .artificial_count(presample, "presample", 0)
  if (!is.character(deterministic) || length(deterministic) != 1 ||
      !deterministic %in% c("none", "const", "trend", "both")) {
    stop("Argument 'deterministic' must be one of 'none', 'const', 'trend' or 'both'.")
  }
  .artificial_flags(list(structural = structural, tvp = tvp, sv = sv, stable = stable))
  .artificial_share(a_zeros, "a_zeros")
  range_a <- .artificial_range(range_a, "range_a")
  range_a0 <- .artificial_range(range_a0, "range_a0")
  range_const <- .artificial_range(range_const, "range_const")
  range_trend <- .artificial_range(range_trend, "range_trend")
  range_variance <- .artificial_range(range_variance, "range_variance", positive = TRUE)
  range_psi <- .artificial_range(range_psi, "range_psi")
  range_variance_state <- .artificial_range(range_variance_state, "range_variance_state", non_negative = TRUE)
  range_variance_sv <- .artificial_range(range_variance_sv, "range_variance_sv", non_negative = TRUE)
  if (!is.numeric(level) || !length(level) %in% c(1, k) || any(!is.finite(level))) {
    stop("Argument 'level' must be a numeric vector with 1 or 'k' finite elements.")
  }
  level <- rep_len(level, k)
  if (structural & any(range_psi != 0)) {
    stop("Argument 'range_psi' must be c(0, 0) for structural models, whose errors are uncorrelated.")
  }

  const <- deterministic %in% c("const", "both")
  if (!const & any(level != 0)) {
    stop("A non-zero 'level' requires an intercept, i.e. 'deterministic' must be 'const' or 'both'.")
  }
  trend <- deterministic %in% c("trend", "both")
  n_lag <- k * p
  m <- n_lag + const + trend
  names_series <- paste0("var", 1:k)

  # Checks whether the reduced form of the coefficients in list 'x' is stable
  is_stable <- function(x) {
    if (!stable | p == 0) {
      return(TRUE)
    }
    lags <- x[["a"]][, 1:n_lag, drop = FALSE]
    if (structural) {
      lags <- forwardsolve(x[["a0"]], lags)
    }
    return(.artificial_var_stable(lags, k, p))
  }

  # Contemporaneous coefficients
  a0 <- NULL
  if (structural) {
    a0 <- diag(1, k)
    if (k > 1) {
      a0[lower.tri(a0)] <- .artificial_draw(k * (k - 1) / 2, range_a0)
    }
    dimnames(a0) <- list(names_series, names_series)
  }

  # Coefficients of the first period
  a <- NULL
  if (m > 0) {
    names_regressors <- c(if (p > 0) paste0(rep(names_series, p), ".", rep(1:p, each = k)),
                          if (const) "const", if (trend) "trend")
    a <- matrix(0, k, m, dimnames = list(names_series, names_regressors))
    if (p > 0) {
      for (i in 1:1000) {
        a[, 1:n_lag] <- .artificial_draw(k * n_lag, range_a, a_zeros)
        if (is_stable(list(a = a, a0 = a0))) {
          break
        }
        if (i == 1000) {
          stop("No stable coefficient matrix was drawn in 1000 attempts. ",
               "Narrow 'range_a', raise 'a_zeros' or set 'stable = FALSE'.")
        }
      }
    }
    if (const) {
      a[, "const"] <- .artificial_draw(k, range_const)
    }
    if (trend) {
      a[, "trend"] <- .artificial_draw(k, range_trend)
    }
  }

  # Error term of the first period
  err <- .artificial_error_init(k, names_series, tvp, sv, range_variance, range_psi,
                                range_variance_state, range_variance_sv)

  # Variances of the state equations
  coef_t <- list(a = a, a0 = a0)
  coef_t <- coef_t[!vapply(coef_t, is.null, logical(1))]
  if (tvp) {
    coef_state <- lapply(coef_t, .artificial_state_variance, range = range_variance_state)
    if (structural) {
      coef_state[["a0"]] <- coef_state[["a0"]] * lower.tri(a0)
    }
  }

  # Initial values of the series at the mean of a stable process without trend
  y <- matrix(0, k, p + presample + nobs)
  if (p > 0 & const & stable) {
    a_reduced <- if (structural) forwardsolve(a0, a) else a
    init_cond <- diag(1, k)
    for (i in 1:p) {
      init_cond <- init_cond - a_reduced[, k * (i - 1) + 1:k, drop = FALSE]
    }
    y[, 1:p] <- solve(init_cond, a_reduced[, n_lag + 1])
  }

  # Generate time series
  hist <- vector("list", nobs)
  for (t in 1:(presample + nobs)) {

    # Position in the returned series
    j <- t - presample

    if (j >= 1) {
      if (tvp & length(coef_t) > 0) {
        coef_t <- .artificial_tvp_step(coef_t, coef_state, is_stable)
      }
      err <- .artificial_error_update(err)
      if (tvp | sv) {
        hist[[j]] <- c(coef_t, err[c("psi_t", "log_omega_t")])
      }
    }

    y_t <- .artificial_error_draw(err)
    if (p > 0) {
      for (i in 1:p) {
        y_t <- y_t + coef_t[["a"]][, k * (i - 1) + 1:k, drop = FALSE] %*% y[, p + t - i]
      }
    }
    if (const) {
      y_t <- y_t + coef_t[["a"]][, "const"]
    }
    if (trend) {
      y_t <- y_t + coef_t[["a"]][, "trend"] * (j - p)
    }
    if (structural) {
      y_t <- forwardsolve(coef_t[["a0"]], y_t)
    }
    y[, p + t] <- y_t
  }

  y <- y[, p + presample + 1:nobs, drop = FALSE] + level

  # The intercept absorbs (A_0 - A_1 - ... - A_p) m of the shift by the level
  shift_const <- function(x) {
    shift <- if (structural) x[["a0"]] %*% level else level
    if (p > 0) {
      for (i in 1:p) {
        shift <- shift - x[["a"]][, k * (i - 1) + 1:k, drop = FALSE] %*% level
      }
    }
    x[["a"]][, "const"] <- x[["a"]][, "const"] + shift
    return(x)
  }
  if (any(level != 0)) {
    a <- shift_const(list(a = a, a0 = a0))[["a"]]
    if (tvp | sv) {
      hist <- lapply(hist, shift_const)
    }
  }

  # Collect true parameters
  errors <- .artificial_error_params(err, hist)
  if (structural) {
    errors[["params"]]["psi_coef"] <- list(NULL)
    errors[["states"]][["psi_state_variance"]] <- NULL
  }

  params <- list(a_coef = if (tvp & m > 0) .artificial_path(hist, "a", a) else a)
  if (structural) {
    params[["a0_coef"]] <- if (tvp) .artificial_path(hist, "a0", a0) else a0
  }
  params <- c(params, errors[["params"]])
  if (tvp) {
    params[["a_state_variance"]] <- coef_state[["a"]]
    params[["a0_state_variance"]] <- coef_state[["a0"]]
  }
  params <- c(params, errors[["states"]])

  data <- stats::ts(t(y))
  dimnames(data) <- list(NULL, names_series)

  result <- list(data = data,
                 params = params)

  return(result)
}
