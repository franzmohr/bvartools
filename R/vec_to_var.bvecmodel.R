#' Transform a VEC Model to a VAR in Levels
#'
#' An object of class \code{'bvecmodel'} is transformed into an object of class
#' \code{'bvarmodel'}, which contains the VAR representation of the model in
#' levels.
#'
#' @param object an object of class \code{'bvecmodel'}, usually, the result of a
#' call to \code{\link{create_bvecmodel}}, optionally already estimated with
#' \code{\link{draw_posterior}}.
#' @param ... arguments passed forward to method.
#'
#' @details A VEC model and its VAR representation in levels are the same model
#' in two parameterisations, so a posterior draw of the one is a posterior draw
#' of the other. The transformation is therefore a change of basis, which is
#' applied draw by draw:
#' \deqn{A_1 = A_0 + \Pi^{y} + \Gamma_1, \quad A_i = \Gamma_i - \Gamma_{i - 1},
#' \quad A_p = -\Gamma_{p - 1}}
#' for the endogenous variables and, analogously, one block later for the
#' unmodelled, non-deterministic variables
#' \deqn{B_0 = \Upsilon_0, \quad B_1 = \Pi^{x} - \Upsilon_0 + \Upsilon_1, \quad
#' B_j = \Upsilon_j - \Upsilon_{j - 1}, \quad B_s = -\Upsilon_{s - 1}.}
#' \eqn{A_0} is the identity matrix unless the model is structural, in which
#' case it is the matrix of contemporaneous coefficients. It enters \eqn{A_1}
#' because it multiplies \eqn{\Delta y_t} in the VEC model and therefore
#' \eqn{y_t} in its VAR representation, which leaves \eqn{A_0 y_{t - 1}} on the
#' right-hand side. The contemporaneous coefficients themselves are not
#' affected by the transformation and are carried over unchanged.
#' Deterministic terms that entered the model unrestricted carry over unchanged.
#' Those restricted to the cointegration space become ordinary regressors of the
#' VAR, entering after the unrestricted ones, and are dropped when the
#' cointegration rank is zero, in which case they had no effect on the model.
#' The coefficient transformation itself is performed by the C++ implementation
#' of the package's model library, which is also used to obtain forecasts of VEC
#' models.
#'
#' The data matrices are reconstructed from the data of the VEC model, so the
#' resulting object covers exactly the same periods. Since a VEC model of lag
#' order \eqn{p} and its VAR representation use the same number of observations,
#' no observations are gained or lost.
#'
#' Note that the covariance matrix of the error term is not affected by the
#' transformation, since both parameterisations describe the same error term.
#'
#' Note also that the elements \code{priors} and \code{initial} of the input
#' object are not carried over, because they describe the coefficients of the
#' VEC model. Use \code{\link{add_priors}} and \code{\link{add_initial_values}}
#' if the resulting model should be estimated in its VAR form. The same applies
#' to the results of variable selection algorithms: a coefficient of the VAR
#' representation is the sum of multiple coefficients of the VEC model, so no
#' single inclusion indicator describes it and the draws of those indicators are
#' dropped.
#'
#' @return An object of class \code{'bvarmodel'}.
#'
#' @examples
#'
#' # Load data
#' data("e6")
#' e6 <- e6 * 100
#'
#' # Generate model
#' model <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # Transform to a VAR model in levels
#' object <- vec_to_var(model)
#'
#' # The result can be used like any other VAR model, e.g. for forecasting
#' object <- add_forecast_input_data(object, n_ahead = 4)
#' object <- add_posterior_forecasts(object)
#'
#' @references
#'
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' @export
#' @method vec_to_var bvecmodel
vec_to_var.bvecmodel <- function(object, ...) {

  specs <- object[["model"]]

  # Dimensions of the VAR representation. Obtained from the model library, so
  # that R and C++ cannot disagree on the shape of the object whose coefficients
  # the same library transforms below.
  var_specs <- .VecToVarSpecification(specs)

  k <- var_specs[["k"]]
  p <- var_specs[["p"]]
  m <- var_specs[["m"]]
  s <- var_specs[["s"]]
  structural <- var_specs[["structural"]]

  rank <- specs[["rank"]]
  if (is.null(rank)) {
    rank <- 0
  }

  # Blocks of the VEC model. 'p' and 's' are lag orders of the VAR
  # representation, so the VEC model has one block fewer of each.
  n_gamma <- max(specs[["p"]] - 1, 0)
  n_upsilon <- ifelse(m > 0, specs[["s"]], 0)
  n_det_ur <- specs[["n"]]
  # Deterministic terms in the cointegration term only affect the model if there
  # is a cointegration term to speak of.
  n_det_r <- ifelse(rank > 0, specs[["n_restricted"]], 0)

  if (m > 0 & n_upsilon == 0) {
    stop("The current values of unmodelled, non-deterministic variables cannot be recovered from a model without differenced regressors.")
  }

  # Data of the VEC model ----

  y_diff <- object[["data"]][["train"]][["y"]]
  w <- object[["data"]][["train"]][["w"]]
  x_diff <- object[["data"]][["train"]][["x"]]

  if (is.null(y_diff) | is.null(w)) {
    stop("Argument 'object' must contain the data of the VEC model in elements 'y' and 'w'.")
  }

  ts_info <- stats::tsp(y_diff)

  y_diff <- as.matrix(y_diff)
  w <- as.matrix(w)
  x_diff <- if (is.null(x_diff)) matrix(NA_real_, nrow(y_diff), 0) else as.matrix(x_diff)

  if (NCOL(w) != k + m + specs[["n_restricted"]]) {
    stop("The error correction term does not have the number of columns implied by the model specification.")
  }
  if (NCOL(x_diff) != k * n_gamma + m * n_upsilon + n_det_ur) {
    stop("The non-cointegration regressors do not have the number of columns implied by the model specification.")
  }

  # Data of the VAR representation ----

  y_names <- specs[["endogen"]]

  # The level of the endogenous variables, recovered from the differences and
  # the levels in the error correction term.
  y <- y_diff + w[, 1:k, drop = FALSE]

  # ... and its lags, obtained by undoing one difference at a time.
  x <- vector("list", p)
  x[[1]] <- w[, 1:k, drop = FALSE]
  if (p > 1) {
    for (i in 2:p) {
      x[[i]] <- x[[i - 1]] - x_diff[, (i - 2) * k + 1:k, drop = FALSE]
    }
  }
  x_names <- paste0(rep(y_names, times = p), ".",
                    rep(.lag_label(1:p, p), each = k))

  # The unmodelled, non-detereministic variables. The VEC model regresses on the
  # current difference, so the levels it implies run from the current period to
  # lag s, with the cointegration term contributing the first lag.
  if (m > 0) {
    exogen_names <- specs[["exogen"]]
    x_exo <- vector("list", s + 1)
    x_exo[[2]] <- w[, k + 1:m, drop = FALSE]
    x_exo[[1]] <- x_exo[[2]] + x_diff[, k * n_gamma + 1:m, drop = FALSE]
    if (s > 1) {
      for (j in 2:s) {
        x_exo[[j + 1]] <- x_exo[[j]] - x_diff[, k * n_gamma + (j - 1) * m + 1:m, drop = FALSE]
      }
    }
    x <- c(x, x_exo)
    x_names <- c(x_names, paste0(rep(exogen_names, times = s + 1), ".l",
                                 rep(.lag_label(0:s, s), each = m)))
  }

  # Deterministic terms. The unrestricted ones keep the order they had, the
  # restricted ones follow, which is the order the transformed coefficients are
  # in as well.
  det_data <- NULL
  if (n_det_ur > 0) {
    det_data <- x_diff[, k * n_gamma + m * n_upsilon + 1:n_det_ur, drop = FALSE]
  }
  if (n_det_r > 0) {
    det_data <- cbind(det_data, w[, k + m + 1:n_det_r, drop = FALSE])
  }
  if (!is.null(det_data)) {
    x <- c(x, list(det_data))
    x_names <- c(x_names, dimnames(det_data)[[2]])
  }

  x <- do.call(cbind, x)

  y <- stats::ts(y, class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- ts_info
  dimnames(y) <- list(NULL, y_names)

  x <- stats::ts(x, class = c("mts", "ts", "matrix"))
  stats::tsp(x) <- ts_info
  dimnames(x) <- list(NULL, x_names)

  if (!is.null(det_data)) {
    det_data <- stats::ts(det_data, class = c("mts", "ts", "matrix"))
    stats::tsp(det_data) <- ts_info
  }

  # SUR form of the regressors, structural data appended as in create_bvarmodel
  z <- kronecker(x, diag(1, k))
  if (structural & k > 1) {
    y_A0 <- kronecker(-y, diag(1, k))
    pos <- NULL
    for (j in 1:k) {
      pos <- c(pos, (j - 1) * k + 1:j)
    }
    z <- cbind(z, y_A0[, -pos])
  }
  dimnames(z) <- NULL

  # Model specification ----

  model <- NULL
  model[["type"]] <- ifelse(structural, "SVAR", ifelse(k == 1, "AR", "VAR"))
  model[["algorithm"]] <- .var_algorithm(specs[["error"]], specs[["tvp"]])
  model[["k"]] <- k
  model[["p"]] <- as.integer(p)
  model[["m"]] <- as.integer(m)
  model[["s"]] <- as.integer(s)
  model[["n"]] <- as.integer(var_specs[["n"]])
  model[["varsel"]] <- var_specs[["varsel"]]
  model[["endogen"]] <- y_names
  if (m > 0) {
    model[["exogen"]] <- specs[["exogen"]]
  }
  model[["structural"]] <- structural
  model[["error"]] <- specs[["error"]]
  model[["tvp"]] <- specs[["tvp"]]
  model[["iterations"]] <- specs[["iterations"]]
  model[["burnin"]] <- specs[["burnin"]]

  result <- list("model" = model,
                 "data" = list("original" = list("endogen" = object[["data"]][["original"]][["endogen"]],
                                                 "exogen" = object[["data"]][["original"]][["exogen"]],
                                                 "deterministic" = det_data),
                               "train" = list("y" = y,
                                              "x" = x,
                                              "z" = z)))

  # Data for forecast simulation is generated in the layout of the VAR
  # representation already, since that is the form a VEC model is forecast in.
  if (!is.null(object[["data"]][["forecast"]])) {
    result[["data"]][["forecast"]] <- object[["data"]][["forecast"]]
    result[["model"]][["h"]] <- specs[["h"]]
  }

  # Posterior draws ----

  if (!is.null(object[["posterior"]])) {

    if (isTRUE(specs[["tvp"]])) {
      stop("Posterior draws of models with time varying parameters cannot be transformed.")
    }

    posterior <- object[["posterior"]]
    mcpar <- attr(posterior[["u_sigma_inv"]][["coeffs"]], "mcpar")

    # The covariance of the error term is shared, not transformed, and so are
    # the log likelihood and any forecast: both are already those of the VAR
    # representation.
    coeffs <- .VecToVarCoefficients(object)[["a"]][["coeffs"]]
    if (is.null(mcpar)) {
      coeffs <- coda::as.mcmc(coeffs)
    } else {
      coeffs <- coda::mcmc(coeffs, start = mcpar[1], end = mcpar[2], thin = mcpar[3])
    }

    posterior[["a"]] <- list("coeffs" = coeffs)
    posterior[["beta"]] <- NULL

    result[["posterior"]] <- posterior
  }

  class(result) <- c("bvarmodel", "list")

  return(result)
}

# Zero-padded lag index, in the format used by create_bvarmodel()
.lag_label <- function(i, max_lag) {
  formatC(as.integer(i), width = max(2, nchar(max_lag)), flag = "0", format = "d")
}

# Name of the posterior simulation algorithm of a VAR model with the given
# specification of the error term, as determined by create_bvarmodel()
.var_algorithm <- function(error, tvp) {

  algorithm <- ifelse(isTRUE(tvp), "Tvp", "Normal")

  if (error == "wishart") {
    algorithm <- paste0(algorithm, "Wishart")
  }
  if (error %in% c("gamma", "gamma+covar")) {
    algorithm <- paste0(algorithm, "Gamma")
  }
  if (error %in% c("sv", "sv+covar")) {
    algorithm <- paste0(algorithm, "Stochvol")
  }

  return(paste0("Var", algorithm))
}
