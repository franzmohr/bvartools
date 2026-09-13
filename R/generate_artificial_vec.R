#' Generate Artificial VEC Data
#'
#' Generates an artificial data set from a vector error correction model with
#' constant or time varying parameters and constant or stochastic volatility,
#' for testing the algorithms of the package.
#'
#' @param nobs number of generated observations in levels. Defaults to 100.
#' @param k number of endogenous variables. Defaults to 3.
#' @param p lag order of the series in the (levels) VAR, so that the VEC model has
#' \eqn{p - 1} lags of differences. Must be at least 1. Defaults to 2.
#' @param r cointegration rank between 0 and \code{k}. Defaults to 1.
#' @param const a character specifying whether a constant term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term is added.
#' @param trend a character specifying whether a linear trend enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no trend is added.
#' @param structural logical specifying whether a structural VEC model with a lower
#' triangular matrix \eqn{A_0} of contemporaneous coefficients is generated.
#' Defaults to \code{FALSE}. See 'Details'.
#' @param tvp logical specifying whether the coefficients, \eqn{A_0} and the
#' elements of \eqn{\Psi} follow random walks. Defaults to \code{FALSE}. See 'Details'.
#' @param sv logical specifying whether the error variances follow a stochastic
#' volatility process. Defaults to \code{FALSE}. See 'Details'.
#' @param range_alpha numeric vector with two elements containing the minimum and
#' maximum value of the loadings \eqn{\alpha}. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_beta numeric vector with two elements containing the minimum and
#' maximum value of the free elements of the cointegration matrix \eqn{\beta}, which
#' belong to endogenous variables. Defaults to \code{c(-1, 1)}.
#' @param range_gamma numeric vector with two elements containing the minimum and
#' maximum value of the coefficients of lagged differences. Defaults to \code{c(-0.5, 0.5)}.
#' @param gamma_zeros numeric between 0 and 1 indicating the share of coefficients
#' of lagged differences, which should be set to zero. Default is \code{0.5}.
#' @param range_a0 numeric vector with two elements containing the minimum and
#' maximum value of the free elements of \eqn{A_0}. Only used if
#' \code{structural = TRUE}. Defaults to \code{c(-0.5, 0.5)}.
#' @param range_const numeric vector with two elements containing the minimum and
#' maximum value of the coefficients of the constant term. Defaults to \code{c(-0.5, 0.5)}.
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
#' @param stable logical specifying whether the coefficients are restricted to a process
#' that is integrated of order one with cointegration rank \code{r}. Defaults to
#' \code{TRUE}. See 'Details'.
#' @param presample numeric specifying the number of observations, which are
#' generated before the first returned observation and then discarded. Defaults to 100.
#'
#' @details The function produces artificial observations for a vector
#' error correction (VEC) model:
#' \deqn{A_{0t} \Delta y_t = \alpha_t \beta_t^{\prime} w_t +
#' \sum_{i=1}^{p-1} \Gamma_{it} \Delta y_{t - i} + C_t d_t + u_t,}
#' where
#' \eqn{\Delta y_t} is a K-dimensional vector of differenced endogenous variables,
#' \eqn{w_t} is the vector of the endogenous variables in levels \eqn{y_{t-1}} and
#' the restricted deterministic terms,
#' \eqn{\alpha_t} is the \eqn{K \times r} matrix of loadings,
#' \eqn{\beta_t} is the cointegration matrix with \eqn{K + N^{R}} rows and \eqn{r} columns,
#' \eqn{\Pi_t = \alpha_t \beta_t^{\prime}},
#' \eqn{\Gamma_{it}} is a \eqn{K \times K} coefficient matrix of lagged differences,
#' \eqn{d_t} is the vector of unrestricted deterministic terms and \eqn{C_t} its
#' coefficient matrix. The error term \eqn{u_t} and \eqn{A_{0t}} are specified as
#' in \code{\link{generate_artificial_var}}: \eqn{u_t \sim N(0, \Sigma_t)} with
#' \eqn{\Psi_t \Sigma_t \Psi_t^{\prime} = \Omega_t}, and \eqn{A_{0t}} is an identity
#' matrix unless \code{structural = TRUE}, where it is lower triangular with ones on
#' its main diagonal and the structural errors are uncorrelated.
#' This is the model that \code{\link{create_bvecmodel}} produces with the same
#' arguments \code{p}, \code{r}, \code{const}, \code{trend} and \code{structural}.
#'
#' The cointegration matrix is normalised so that its first \eqn{r} rows are an identity
#' matrix. Its other rows for endogenous variables are drawn from \code{range_beta}, and
#' its rows for a restricted constant or trend from \code{range_const} and
#' \code{range_trend}. Since \code{\link{add_posterior_coefficients}} may use another
#' normalisation, estimates are best compared with the true parameters through
#' \eqn{\Pi_t}, which does not depend on it.
#'
#' The linear trend takes the value 1 in period \eqn{p + 1} of the returned series,
#' which is the first period that \code{\link{create_bvecmodel}} uses for estimation
#' with lag order \eqn{p}.
#'
#' If \code{tvp = TRUE}, the loadings, the free elements of \eqn{\beta_t}, the
#' coefficients of lagged differences and unrestricted deterministic terms, the free
#' elements of \eqn{A_{0t}} and the free elements of \eqn{\Psi_t} follow random walks,
#' whose innovation variances are drawn from \code{range_variance_state}. Coefficients
#' that are zero in the first period stay zero. The random walk of the normalised
#' \eqn{\beta_t} differs from the state equation of the unnormalised cointegration
#' matrix in Koop et al. (2011), which \code{\link{create_bvecmodel}} uses for
#' estimation. If \code{sv = TRUE}, the log-volatilities follow random walks as in
#' \code{\link{generate_artificial_var}}. During the presample the parameters are held
#' at their values of the first period.
#'
#' If \code{stable = TRUE}, the coefficients are drawn until the VAR representation of
#' the reduced form of the process has exactly \eqn{K - r} unit roots and all other
#' roots inside the unit circle, so that the series are integrated of order one and
#' \eqn{\beta_t^{\prime} w_t} is stationary. With \code{r = 0} the series are
#' integrated but not cointegrated, with \code{r = k} they are stationary. For
#' \code{tvp = TRUE} this holds in every period: an innovation, which violates it, is
#' drawn again, and if no admissible innovation is found in 100 attempts, the
#' coefficients keep the values of the previous period.
#'
#' @references
#' Johansen, S. (1995). \emph{Likelihood-based inference in cointegrated vector
#' autoregressive models}. Oxford: Oxford University Press.
#'
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#'
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' @returns A list with the elements
#' \describe{
#'   \item{\code{data}}{a \eqn{T \times K} time-series object of the artificial
#'   series in levels, named \code{var1}, \code{var2} etc.}
#'   \item{\code{params}}{a list of the true parameters with the elements
#'   \describe{
#'     \item{\code{alpha}}{the \eqn{K \times r} matrix of loadings \eqn{\alpha}.}
#'     \item{\code{beta}}{the \eqn{(K + N^{R}) \times r} cointegration matrix \eqn{\beta}.}
#'     \item{\code{pi}}{the \eqn{K \times (K + N^{R})} matrix \eqn{\Pi = \alpha \beta^{\prime}},
#'     with the columns in the order of the cointegration variables of
#'     \code{\link{create_bvecmodel}}.}
#'     \item{\code{gamma}}{the \eqn{K \times K(p - 1)} matrix \eqn{(\Gamma_1, \dots, \Gamma_{p-1})}.}
#'     \item{\code{c}}{the \eqn{K \times N^{UR}} matrix \eqn{C} of unrestricted deterministic terms.}
#'     \item{\code{a0_coef}}{for \code{structural = TRUE}, the \eqn{K \times K} matrix \eqn{A_0}.}
#'     \item{\code{psi_coef}, \code{u_omega}, \code{u_sigma}}{the parameters of the error
#'     term as in \code{\link{generate_artificial_var}}.}
#'     \item{\code{alpha_state_variance}, \code{beta_state_variance},
#'     \code{gamma_state_variance}, \code{c_state_variance}, \code{a0_state_variance},
#'     \code{psi_state_variance}}{for \code{tvp = TRUE}, the variances of the state
#'     equations of the respective parameters.}
#'     \item{\code{u_state_variance}}{for \code{sv = TRUE}, the \eqn{K} variances of
#'     the state equations of the log-volatilities.}
#'   }}
#' }
#' Parameters that are not part of the model are \code{NULL}. Time varying parameters
#' are \eqn{\dots \times T} arrays with one matrix per period.
#'
#' @examples
#'
#' # Set seed of RNG
#' set.seed(1)
#'
#' # Three series with one cointegration relationship
#' dt <- generate_artificial_vec(nobs = 200, k = 3, p = 2, r = 1)
#' dt[["params"]][["pi"]]
#'
#' # Restricted constant and time varying cointegration with stochastic volatility
#' dt <- generate_artificial_vec(nobs = 200, k = 2, p = 1, r = 1, const = "restricted",
#'                               tvp = TRUE, sv = TRUE, range_variance_state = c(0.001, 0.001))
#'
#' # Path of the cointegration coefficient of the second variable
#' plot(dt[["params"]][["beta"]][2, 1, ], type = "l")
#'
#' # Structural model with an unrestricted constant
#' dt <- generate_artificial_vec(nobs = 200, k = 3, p = 2, r = 1, const = "unrestricted",
#'                               structural = TRUE)
#'
#' @family artificial data
#' @export
generate_artificial_vec <- function(nobs = 100, k = 3, p = 2, r = 1,
                                    const = NULL,
                                    trend = NULL,
                                    structural = FALSE,
                                    tvp = FALSE,
                                    sv = FALSE,
                                    range_alpha = c(-0.5, 0.5),
                                    range_beta = c(-1, 1),
                                    range_gamma = c(-0.5, 0.5),
                                    gamma_zeros = 0.5,
                                    range_a0 = c(-0.5, 0.5),
                                    range_const = c(-0.5, 0.5),
                                    range_trend = c(-0.01, 0.01),
                                    range_variance = c(1, 1),
                                    range_psi = c(0, 0),
                                    range_variance_state = c(0.0001, 0.0001),
                                    range_variance_sv = c(0.01, 0.01),
                                    stable = TRUE,
                                    presample = 100) {

  # Basic checks
  .artificial_count(nobs, "nobs", 1)
  .artificial_count(k, "k", 1)
  .artificial_count(p, "p", 1)
  .artificial_count(r, "r", 0)
  .artificial_count(presample, "presample", 0)
  if (r > k) {
    stop("Argument 'r' must not be larger than 'k'.")
  }
  for (i in c("const", "trend")) {
    x <- get(i)
    if (!is.null(x) && (!is.character(x) || length(x) != 1 || !x %in% c("restricted", "unrestricted"))) {
      stop(paste0("Argument '", i, "' must be NULL, 'restricted' or 'unrestricted'."))
    }
  }
  .artificial_flags(list(structural = structural, tvp = tvp, sv = sv, stable = stable))
  .artificial_share(gamma_zeros, "gamma_zeros")
  range_alpha <- .artificial_range(range_alpha, "range_alpha")
  range_beta <- .artificial_range(range_beta, "range_beta")
  range_gamma <- .artificial_range(range_gamma, "range_gamma")
  range_a0 <- .artificial_range(range_a0, "range_a0")
  range_const <- .artificial_range(range_const, "range_const")
  range_trend <- .artificial_range(range_trend, "range_trend")
  range_variance <- .artificial_range(range_variance, "range_variance", positive = TRUE)
  range_psi <- .artificial_range(range_psi, "range_psi")
  range_variance_state <- .artificial_range(range_variance_state, "range_variance_state", non_negative = TRUE)
  range_variance_sv <- .artificial_range(range_variance_sv, "range_variance_sv", non_negative = TRUE)
  if (structural & any(range_psi != 0)) {
    stop("Argument 'range_psi' must be c(0, 0) for structural models, whose errors are uncorrelated.")
  }

  const_r <- identical(const, "restricted")
  const_ur <- identical(const, "unrestricted")
  trend_r <- identical(trend, "restricted")
  trend_ur <- identical(trend, "unrestricted")
  names_det_r <- c(if (const_r) "const", if (trend_r) "trend")
  names_det_ur <- c(if (const_ur) "const", if (trend_ur) "trend")
  if (r == 0 & length(names_det_r) > 0) {
    stop("Restricted deterministic terms require a cointegration rank 'r' of at least 1.")
  }
  names_series <- paste0("var", 1:k)

  # Checks whether the reduced form of the coefficients in list 'x' is integrated
  # of order one with cointegration rank r
  is_stable <- function(x) {
    if (!stable) {
      return(TRUE)
    }
    alpha <- x[["alpha"]]
    gamma <- x[["gamma"]]
    if (structural) {
      if (r > 0) {
        alpha <- forwardsolve(x[["a0"]], alpha)
      }
      if (p > 1) {
        gamma <- forwardsolve(x[["a0"]], gamma)
      }
    }
    return(.artificial_vec_stable(alpha, x[["beta"]], gamma, k, p, r))
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
  alpha <- NULL
  beta <- NULL
  gamma <- NULL
  c_det <- NULL
  for (i in 1:1000) {
    if (r > 0) {
      alpha <- matrix(.artificial_draw(k * r, range_alpha), k, r)
      beta <- rbind(diag(1, r), matrix(.artificial_draw((k - r) * r, range_beta), k - r, r))
    }
    if (p > 1) {
      gamma <- matrix(.artificial_draw(k * k * (p - 1), range_gamma, gamma_zeros), k)
    }
    if (is_stable(list(alpha = alpha, beta = beta, gamma = gamma, a0 = a0))) {
      break
    }
    if (i == 1000) {
      stop("No coefficients of a process with cointegration rank ", r, " were drawn in 1000 attempts. ",
           "Narrow 'range_alpha', 'range_beta' or 'range_gamma', or set 'stable = FALSE'.")
    }
  }
  if (r > 0) {
    if (length(names_det_r) > 0) {
      beta_det <- c(if (const_r) .artificial_draw(r, range_const),
                    if (trend_r) .artificial_draw(r, range_trend))
      beta <- rbind(beta, matrix(beta_det, ncol = r, byrow = TRUE))
    }
    names_ect <- paste0("ect", 1:r)
    dimnames(alpha) <- list(names_series, names_ect)
    dimnames(beta) <- list(c(paste0("l.", names_series), names_det_r), names_ect)
  }
  if (p > 1) {
    lags <- formatC(rep(1:(p - 1), each = k), width = max(2, nchar(p - 1)), flag = "0")
    dimnames(gamma) <- list(names_series, paste0("d.", rep(names_series, p - 1), ".l", lags))
  }
  if (length(names_det_ur) > 0) {
    c_det <- matrix(c(if (const_ur) .artificial_draw(k, range_const),
                      if (trend_ur) .artificial_draw(k, range_trend)),
                    k, dimnames = list(names_series, names_det_ur))
  }

  # Error term of the first period
  err <- .artificial_error_init(k, names_series, tvp, sv, range_variance, range_psi,
                                range_variance_state, range_variance_sv)

  # Variances of the state equations. The normalised rows of beta stay fixed.
  coef_t <- list(alpha = alpha, beta = beta, gamma = gamma, c = c_det, a0 = a0)
  coef_t <- coef_t[!vapply(coef_t, is.null, logical(1))]
  if (tvp) {
    coef_state <- lapply(coef_t, .artificial_state_variance, range = range_variance_state)
    if (r > 0) {
      coef_state[["beta"]][1:r, ] <- 0
    }
    if (structural) {
      coef_state[["a0"]] <- coef_state[["a0"]] * lower.tri(a0)
    }
  }

  # Generate time series
  y <- matrix(0, k, p + presample + nobs)
  hist <- vector("list", nobs)
  for (t in 1:(presample + nobs)) {

    # Position in the returned series and value of the trend
    j <- t - presample
    trend_t <- j - p

    if (j >= 1) {
      if (tvp & length(coef_t) > 0) {
        coef_t <- .artificial_tvp_step(coef_t, coef_state, is_stable)
      }
      err <- .artificial_error_update(err)
      if (tvp | sv) {
        hist[[j]] <- c(coef_t, err[c("psi_t", "log_omega_t")])
      }
    }

    dy <- .artificial_error_draw(err)
    if (r > 0) {
      w <- c(y[, p + t - 1], if (const_r) 1, if (trend_r) trend_t)
      dy <- dy + coef_t[["alpha"]] %*% crossprod(coef_t[["beta"]], w)
    }
    if (p > 1) {
      for (i in 1:(p - 1)) {
        dy <- dy + coef_t[["gamma"]][, (i - 1) * k + 1:k, drop = FALSE] %*% (y[, p + t - i] - y[, p + t - i - 1])
      }
    }
    if (const_ur) {
      dy <- dy + coef_t[["c"]][, "const"]
    }
    if (trend_ur) {
      dy <- dy + coef_t[["c"]][, "trend"] * trend_t
    }
    if (structural) {
      dy <- forwardsolve(coef_t[["a0"]], dy)
    }
    y[, p + t] <- y[, p + t - 1] + dy
  }

  y <- y[, p + presample + 1:nobs, drop = FALSE]

  # Collect true parameters
  errors <- .artificial_error_params(err, hist)
  if (structural) {
    errors[["params"]]["psi_coef"] <- list(NULL)
    errors[["states"]][["psi_state_variance"]] <- NULL
  }

  path <- function(name, x) {
    if (tvp & !is.null(x)) .artificial_path(hist, name, x) else x
  }

  pi <- NULL
  if (r > 0) {
    if (tvp) {
      pi <- array(NA_real_, c(k, nrow(beta), nobs), dimnames = list(names_series, rownames(beta), NULL))
      for (j in 1:nobs) {
        pi[, , j] <- hist[[j]][["alpha"]] %*% t(hist[[j]][["beta"]])
      }
    } else {
      pi <- alpha %*% t(beta)
    }
  }

  params <- list(alpha = path("alpha", alpha),
                 beta = path("beta", beta),
                 pi = pi,
                 gamma = path("gamma", gamma),
                 c = path("c", c_det))
  if (structural) {
    params[["a0_coef"]] <- path("a0", a0)
  }
  params <- c(params, errors[["params"]])
  if (tvp) {
    for (i in names(coef_state)) {
      params[[paste0(i, "_state_variance")]] <- coef_state[[i]]
    }
  }
  params <- c(params, errors[["states"]])

  data <- stats::ts(t(y))
  dimnames(data) <- list(NULL, names_series)

  result <- list(data = data,
                 params = params)

  return(result)
}
