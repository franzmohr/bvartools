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
#' the posterior draws of \eqn{alpha} and \eqn{beta} to the original scale
#' of the input data, based on the scaling factors that are stored in attribute
#' \code{"scale"} of \code{object$data$train$w}.
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
  
  if (r > 0) {
    draws <- nrow(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    k <- object[["model"]][["k"]]
    n_alpha <- k * r
    k_ect <- ncol(object[["data"]][["train"]][["w"]])
    n_beta <- k_ect * r
    pos_alpha <- 1:n_alpha
    
    for (i in 1:draws) {
      alpha_i <- rescale_matrix %*% matrix(object[["posterior"]][["a"]][["coeffs"]][i, pos_alpha], k)
      object[["posterior"]][["a"]][["coeffs"]][draw, pos_alpha] <- alpha_i
      
      beta_i <- rescale_matrix_inv %*% matrix(object[["posterior"]][["beta"]][["coeffs"]][i,], k_ect)
      object[["posterior"]][["beta"]][["coeffs"]][draw,] <- beta_i
    }
  }
  
  
  return(object)
}