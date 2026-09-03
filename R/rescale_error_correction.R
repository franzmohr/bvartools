#' Rescale Error Correction
#' 
#' Rescales the series in the error correction series.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
rescale_error_correction <- function (object, ...) {
  UseMethod("rescale_error_correction")
}



#' Rescale Error Correction
#' 
#' Rescales the series in the error correction series in an object of class
#' 'bvecmodel'.
#' 
#' @param object object of class 'bvecmodel'.
#' @param ... arguments passed forward to method.
#' 
#' @details The function transforms element \code{object$data$train$w} and
#' the posterior draws of \eqn{\beta} to the original scale of the input data,
#' based on the scaling factors that are stored in attribute \code{"scale"} of
#' \code{object$data$train$w}.
#' 
#' If \eqn{D} is the diagonal matrix of scaling factors, function
#' \code{\link{scale_error_correction}} replaced \eqn{w_t} by
#' \eqn{D^{-1} w_t}, so that the estimated error correction term is
#' \eqn{\alpha \beta^{\prime} D^{-1} w_t}. The draws of \eqn{\beta} are
#' therefore multiplied by \eqn{D^{-1}}, while the draws of \eqn{\alpha} are
#' not affected by the transformation and are carried over unchanged.
#' 
#' @return An object of class 'bvecmodel'.
#' 
#' @export
#' @method rescale_error_correction bvecmodel
rescale_error_correction.bvecmodel <- function(object, ...) {
  
  # Get rescale factor
  if (!is.null(attr(object[["data"]][["train"]][["w"]], "scale"))) {
    rescale_factor <- attr(object[["data"]][["train"]][["w"]], "scale")
    rescale_matrix <- diag(rescale_factor)
    rescale_matrix_inv <- diag(1 / rescale_factor)
  } else {
    stop("Element 'object$data$train$w' does not have an attribute 'scale'.")
  }
  
  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Model does not seem to contain posterior draws.")
  }
  
  # Input data
  object[["data"]][["train"]][["w"]][] <- t(rescale_matrix %*% t(object[["data"]][["train"]][["w"]]))
  
  r <- object[["model"]][["rank"]]
  
  # Only the draws of beta are affected. Since the error correction term of the
  # estimated model is alpha %*% t(beta) %*% solve(D) %*% w, the coefficients on
  # the original scale are obtained by multiplying beta by solve(D), which
  # leaves alpha unchanged.
  if (r > 0) {
    draws <- nrow(object[["posterior"]][["beta"]][["coeffs"]])
    k_ect <- ncol(object[["data"]][["train"]][["w"]])
    
    for (draw in 1:draws) {
      beta_i <- rescale_matrix_inv %*% matrix(object[["posterior"]][["beta"]][["coeffs"]][draw,], k_ect)
      object[["posterior"]][["beta"]][["coeffs"]][draw,] <- beta_i
    }
  }
  
  # The data are back on their original scale, so there is nothing left to
  # rescale. Dropping the attribute makes a second call fail instead of
  # silently transforming the model a second time.
  attr(object[["data"]][["train"]][["w"]], "scale") <- NULL
  
  return(object)
}