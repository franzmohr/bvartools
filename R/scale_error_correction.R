#' Scale Error Correction
#'
#' Scales the series in the error correction series.
#'
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{scale_error_correction.bvecmodel}}.
#'
#' @export
scale_error_correction <- function (object, ...) {
  UseMethod("scale_error_correction")
}



#' Scale Error Correction
#'
#' Scales, centres, or centres and scales the series in the error correction
#' term of an object of class 'bvecmodel'.
#'
#' @param object object of class 'bvecmodel'.
#' @param scale logical. If \code{TRUE}, the default, the series are divided by
#' scaling factors. See 'Details'.
#' @param centre logical. If \code{TRUE}, the sample mean of each stochastic
#' series is subtracted from it. Default is \code{FALSE}. See 'Details'.
#' @param ... arguments passed forward to method.
#'
#' @details The function transforms element \code{object$data$train$w}, and
#' \code{\link{rescale_error_correction}} undoes the transformation once the
#' posterior draws are in.
#'
#' With \code{scale = TRUE}, stochastic series are divided by the standard
#' deviation of the corresponding differenced series. If the time-series object
#' contains a column named \code{"trend"}, this series is divided by its own
#' standard deviation, i.e. in levels. The scaling factors are stored as a new
#' attribute of \code{object$data$train$w} named \code{"scale"}.
#'
#' With \code{centre = TRUE}, the sample mean of each stochastic series is
#' subtracted from it. The deterministic terms restricted to the cointegration
#' space, which are the last columns of the term, are left as they are. The
#' means are stored as attribute \code{"centre"}, with zeros for the
#' deterministic terms. If both arguments are \code{TRUE}, the series are
#' centred first and scaled second, so that \eqn{w_t} becomes
#' \eqn{D^{-1} (w_t - m)}.
#'
#' Centring requires an unrestricted constant, which takes up what the series
#' lose: \eqn{\Pi w_t + c = \Pi (w_t - m) + (c + \Pi m)}. For a model with
#' constant coefficients it is therefore a reparameterisation that leaves the
#' likelihood and the priors on \eqn{\alpha} and \eqn{\beta} as they are and
#' changes only what the prior of the constant refers to. For a model with
#' time varying cointegration vectors it changes more. A step \eqn{\eta_t} of
#' their state equation moves \eqn{\beta_t^{\prime} w_t} by
#' \eqn{\eta_t^{\prime} w_t}, which for series far from zero acts as a random
#' walk intercept that no prior on the deterministic terms controls; see
#' section 'Prior on the cointegration space' of \code{\link{cointspace_prior}}.
#' On centred series the same step moves the term by
#' \eqn{\eta_t^{\prime} (w_t - m)}, and the intercept is left to the constant
#' and its own prior. That holds for the estimation. The forecast of such a
#' model, which needs \code{\link{rescale_error_correction}} first, can only hold
#' the cointegration vectors at their last value (\code{forecast_states = "hold"}
#' in \code{\link{add_posterior_forecasts.bvecmodel}}), since after rescaling a
#' step of them acts on the series as they are.
#'
#' Starting values of the constant in \code{object$initial} are shifted by
#' \eqn{\Pi m}, computed from the starting values of the loadings and of
#' \eqn{\beta}, so that the chain starts from the model it simulates.
#'
#' Neither transformation can be applied twice, and neither can be added to a
#' model that already carries the other: \code{rescale_error_correction} has
#' to be called in between.
#'
#' @return An object of class 'bvecmodel'.
#'
#' @export
#' @method scale_error_correction bvecmodel
scale_error_correction.bvecmodel <- function(object, scale = TRUE, centre = FALSE, ...) {

  is_flag <- function(x) {is.logical(x) && length(x) == 1 && !is.na(x)}
  if (!is_flag(scale) || !is_flag(centre)) {
    stop("Arguments 'scale' and 'centre' must be TRUE or FALSE.")
  }
  if (!scale && !centre) {
    stop("At least one of the arguments 'scale' and 'centre' must be TRUE.")
  }

  w <- object[["data"]][["train"]][["w"]]

  # Scaling twice is a no-op on the numbers -- dividing by the standard
  # deviation of the differences makes that standard deviation one, so the
  # second pass divides by one -- but it recomputes the factors from the series
  # it is given and overwrites the stored ones with those ones. The first
  # scaling would then be irreversible. Centring twice is the same: the second
  # pass finds means of zero and overwrites the ones that are needed on the way
  # back. Refused, the way a second rescaling is by the attributes being dropped.
  if (!is.null(attr(w, "scale")) || !is.null(attr(w, "centre"))) {
    stop("The series in the error correction term are already scaled or ",
         "centred. Use 'rescale_error_correction' to put them back on the scale ",
         "of the input data first.")
  }

  # A cointegration space prior other than a multiple of the identity has a
  # direction, and that direction is stated in the units of 'w' it was built
  # against -- for coint$p_tau_i = "ml" the maximum likelihood estimate on the
  # series as they were. Scaling afterwards would leave it pointing somewhere
  # else without any sign of it. That holds for the constant model's p_tau_inv
  # and for the time varying model's transition p_tau alike. Centring does not
  # move it: the estimate is computed after the unrestricted constant has been
  # partialled out, which takes the means with it.
  if (scale) {
    has_direction <- function(m) {
      !is.null(m) && !isTRUE(all.equal(m, diag(m[1, 1], nrow(m)), check.attributes = FALSE))
    }
    # The shrinkage of a constant model's prior is not free of the scale
    # either, whatever P_tau is: it sets the prior of alpha given beta, and
    # rescaling 'w' rescales beta, so the same v_inv then means a different
    # prior on alpha beta'. For coint$v_i = "ml" it was also computed from the
    # series as they were. Only v_inv = 0, the uniform prior, is untouched.
    beta_prior <- object[["priors"]][["beta"]]
    informative <- is.null(beta_prior[["rho"]]) &&
      length(beta_prior[["v_inv"]]) == 1 && isTRUE(beta_prior[["v_inv"]] > 0)
    if (informative || has_direction(beta_prior[["p_tau_inv"]]) ||
        has_direction(beta_prior[["p_tau"]])) {
      stop("The model already has a cointegration space prior that depends on the ",
           "scale of the error correction term. Call 'scale_error_correction' ",
           "before 'add_priors'.")
    }
  }

  tt <- nrow(w)
  values <- matrix(as.numeric(w), tt)

  if (centre) {

    if (is.null(.ect_constant_positions(object))) {
      stop("Centring the error correction term requires an unrestricted constant, ",
           "which takes up the means of the series. Use const = \"unrestricted\" ",
           "or leave 'centre' FALSE.")
    }

    # The deterministic terms restricted to the cointegration space come last.
    # A restricted constant would be a column of zeros once centred.
    n_restricted <- object[["model"]][["n_restricted"]]
    if (is.null(n_restricted)) {
      n_restricted <- 0
    }
    stochastic <- seq_len(ncol(w) - n_restricted)

    means <- rep(0, ncol(w))
    means[stochastic] <- colMeans(values[, stochastic, drop = FALSE])
    names(means) <- dimnames(w)[[2]]

    # The starting values, while beta is still in the units of the series as
    # they are.
    object <- .shift_initial_ect_constant(object, means)

    values <- values - rep(means, each = tt)
    attr(object[["data"]][["train"]][["w"]], "centre") <- means
  }

  if (scale) {

    rescale_factors <- rep(1, ncol(w))
    pos_non_trend <- which(!dimnames(w)[[2]] %in% c("const", "trend"))
    rescale_factors[pos_non_trend] <- apply(diff(w[, pos_non_trend, drop = FALSE]), 2, sd)

    if ("trend" %in% dimnames(w)[[2]]) {
      pos_trend <- which(dimnames(w)[[2]] == "trend")
      rescale_factors[pos_trend] <- sd(w[, "trend"])
    }

    names(rescale_factors) <- dimnames(w)[[2]]
    attr(object[["data"]][["train"]][["w"]], "scale") <- rescale_factors

    values <- values / rep(rescale_factors, each = tt)
  }

  object[["data"]][["train"]][["w"]][] <- values

  return(object)
}



# Positions of the unrestricted constant among the coefficients of one period.
#
# The coefficients of a period are vec(alpha, C) with the loadings first and the
# columns of the non-cointegration regressors 'x' after them, so the k elements
# of the constant are the column of "const" in 'x', counted after the k * r
# loadings. NULL if the model has no unrestricted constant.
.ect_constant_positions <- function(object) {

  x <- object[["data"]][["train"]][["x"]]
  if (is.null(x)) {
    return(NULL)
  }
  column <- which(dimnames(x)[[2]] == "const")
  if (length(column) != 1) {
    return(NULL)
  }

  k <- object[["model"]][["k"]]
  r <- object[["model"]][["rank"]]
  if (is.null(r)) {
    r <- 0
  }

  k * r + (column - 1) * k + seq_len(k)
}



# Add alpha beta' location to the unrestricted constant of every row of 'a'.
#
# 'a' holds the coefficients of one or more periods per row -- draws, or a
# single row of starting values -- and 'beta' the cointegration vectors of the
# same periods. A model whose coefficients vary over time has the loadings and
# the cointegration vectors of each period, so the shift is formed period by
# period. 'location' is in the units of the series the beta in 'beta' belongs
# to. The number of coefficients of a period is that of 'data$train$z'.
.shift_ect_constant <- function(object, a, beta, location) {

  k <- object[["model"]][["k"]]
  r <- object[["model"]][["rank"]]
  k_ect <- length(location)
  n_coef <- NCOL(object[["data"]][["train"]][["z"]])
  constant <- .ect_constant_positions(object)

  if (is.null(r) || r == 0 || is.null(beta)) {
    return(a)
  }
  if (n_coef == 0) {
    n_coef <- k * r + k * NCOL(object[["data"]][["train"]][["x"]])
  }

  a_periods <- ncol(a) / n_coef
  beta_periods <- ncol(beta) / (k_ect * r)
  if (a_periods != round(a_periods) || beta_periods != round(beta_periods) ||
      (beta_periods != 1 && beta_periods != a_periods)) {
    stop("The coefficients and the cointegration vectors do not have a matching number of periods.")
  }

  for (period in seq_len(a_periods)) {
    offset_a <- (period - 1) * n_coef
    offset_beta <- (min(period, beta_periods) - 1) * k_ect * r

    shift <- matrix(0, nrow(a), k)
    for (l in seq_len(r)) {
      beta_l <- beta[, offset_beta + (l - 1) * k_ect + seq_len(k_ect), drop = FALSE]
      alpha_l <- a[, offset_a + (l - 1) * k + seq_len(k), drop = FALSE]
      shift <- shift + alpha_l * as.numeric(beta_l %*% location)
    }

    a[, offset_a + constant] <- a[, offset_a + constant] + shift
  }

  a
}



# Shift the starting values of the constant -- the path and, for a model with
# time varying parameters, the state before the sample -- by alpha beta'
# location, with the starting values of alpha and beta.
.shift_initial_ect_constant <- function(object, location) {

  initial <- object[["initial"]]
  if (is.null(initial[["a"]]) || is.null(initial[["beta"]])) {
    return(object)
  }

  pairs <- list(c("a", "beta"), c("a_init", "beta_init"))
  for (pair in pairs) {
    a <- initial[[pair[1]]]
    beta <- initial[[pair[2]]]
    if (is.null(a) || is.null(beta)) {
      next
    }
    shifted <- .shift_ect_constant(object, matrix(as.numeric(a), 1),
                                   matrix(as.numeric(beta), 1), location)
    a[] <- as.numeric(shifted)
    object[["initial"]][[pair[1]]] <- a
  }

  return(object)
}
