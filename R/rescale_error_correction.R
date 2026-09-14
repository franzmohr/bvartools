#' Rescale Error Correction
#'
#' Rescales the series in the error correction series.
#'
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{rescale_error_correction.bvecmodel}}.
#'
#' @export
rescale_error_correction <- function (object, ...) {
  UseMethod("rescale_error_correction")
}



#' Rescale Error Correction
#'
#' Puts the series in the error correction term of an object of class
#' 'bvecmodel', and the posterior draws that belong to them, back in terms of
#' the input data.
#'
#' @param object object of class 'bvecmodel'.
#' @param ... arguments passed forward to method.
#'
#' @details The function transforms element \code{object$data$train$w}, the
#' posterior draws and the starting values back to the series of the input
#' data, based on the attributes \code{"scale"} and \code{"centre"} that
#' \code{\link{scale_error_correction}} stored in \code{object$data$train$w},
#' and drops the attributes.
#'
#' If \eqn{D} is the diagonal matrix of scaling factors, function
#' \code{\link{scale_error_correction}} replaced \eqn{w_t} by
#' \eqn{D^{-1} w_t}, so that the estimated error correction term is
#' \eqn{\alpha \beta^{\prime} D^{-1} w_t}. The draws of \eqn{\beta} are
#' therefore multiplied by \eqn{D^{-1}}, while the draws of \eqn{\alpha} are
#' not affected by the transformation and are carried over unchanged.
#'
#' If the series were centred on their means \eqn{m}, the estimated term is
#' \eqn{\alpha \beta^{\prime} D^{-1} (w_t - m)}, with \eqn{D} the identity
#' matrix if they were not scaled. The draws of the unrestricted constant are
#' therefore shifted by \eqn{-\alpha \beta^{\prime} D^{-1} m}, draw by draw
#' and, for a model with time varying parameters, period by period with the
#' loadings and the cointegration vectors of that period. The fitted values and
#' the log-likelihood of every draw are unchanged by that. The starting values
#' of the constant are shifted in the same way. If variable selection covered
#' the deterministic terms, the draws of the inclusion indicators of the
#' constant are carried over as they are, so they describe the constant of the
#' centred series.
#'
#' @return An object of class 'bvecmodel'.
#'
#' @export
#' @method rescale_error_correction bvecmodel
rescale_error_correction.bvecmodel <- function(object, ...) {

  w <- object[["data"]][["train"]][["w"]]
  rescale_factor <- attr(w, "scale")
  means <- attr(w, "centre")

  if (is.null(rescale_factor) && is.null(means)) {
    stop("Element 'object$data$train$w' does not have an attribute 'scale' or 'centre'.")
  }

  k_ect <- ncol(w)
  tt <- nrow(w)
  factors <- if (is.null(rescale_factor)) rep(1, k_ect) else rescale_factor
  # Given an explicit size, because for a single factor diag() would read its
  # argument as the size of an identity matrix rather than as the diagonal to
  # build.
  rescale_matrix_inv <- diag(1 / factors, nrow = k_ect)

  # Draws are not required. Putting the series back on the scale of the data is
  # the part that always applies, and a model can be scaled before it has been
  # estimated -- that is the state one is exported in for an external sampler,
  # and the state a sub-model is left in when its run did not produce draws.
  # Only the transformation of the coefficients below needs them.
  has_draws <- !is.null(object[["posterior"]][["beta"]][["coeffs"]])

  r <- object[["model"]][["rank"]]
  if (is.null(r)) {
    r <- 0
  }

  # The constant first, while beta is still in the units of the series the
  # sampler saw: the estimated term is alpha beta' D^-1 (w - m), so the constant
  # of the series as they are is the estimated one minus alpha beta' D^-1 m.
  if (!is.null(means) && r > 0) {
    if (is.null(.ect_constant_positions(object))) {
      stop("The error correction term is centred, but the model has no unrestricted ",
           "constant to take the means back from.")
    }
    location <- -means / factors
    if (has_draws && !is.null(object[["posterior"]][["a"]][["coeffs"]])) {
      object[["posterior"]][["a"]][["coeffs"]] <-
        .shift_ect_constant(object, object[["posterior"]][["a"]][["coeffs"]],
                            object[["posterior"]][["beta"]][["coeffs"]], location)
    }
    object <- .shift_initial_ect_constant(object, location)
  }

  # Input data
  values <- matrix(as.numeric(w), tt) * rep(factors, each = tt)
  if (!is.null(means)) {
    values <- values + rep(means, each = tt)
  }
  object[["data"]][["train"]][["w"]][] <- values

  # Only the draws of beta are scaled. Since the error correction term of the
  # estimated model is alpha %*% t(beta) %*% solve(D) %*% w, the coefficients on
  # the original scale are obtained by multiplying beta by solve(D), which
  # leaves alpha unchanged.
  if (!is.null(rescale_factor) && r > 0 && has_draws) {
    draws <- nrow(object[["posterior"]][["beta"]][["coeffs"]])

    for (draw in 1:draws) {
      beta_i <- rescale_matrix_inv %*% matrix(object[["posterior"]][["beta"]][["coeffs"]][draw,], k_ect)
      object[["posterior"]][["beta"]][["coeffs"]][draw,] <- beta_i
    }
  }

  # The data are back on their original scale, so there is nothing left to
  # rescale. Dropping the attributes makes a second call fail instead of
  # silently transforming the model a second time.
  attr(object[["data"]][["train"]][["w"]], "scale") <- NULL
  attr(object[["data"]][["train"]][["w"]], "centre") <- NULL

  return(object)
}
