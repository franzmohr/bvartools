#' Scale Error Correction
#' 
#' Scales the series in the error correction series.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @export
scale_error_correction <- function (object, ...) {
  UseMethod("scale_error_correction")
}



#' Scale Error Correction
#' 
#' Scales the series in the error correction series in an object of class
#' 'bvecmodel'.
#' 
#' @param object object of class 'bvecmodel'.
#' @param ... arguments passed forward to method.
#' 
#' @details The function transforms element \code{object$data$train$w}. For
#' stochastic variables, the series are divided by the standard deviation of
#' the corresponding differenced series. If the time-series object contains
#' a column named \code{"trend"}, this series is divided by its own standard
#' deviation, i.e. in levels. The scaling factors are stored as a new attribute
#' \code{object$data$train$w} named \code{"scale"}. 
#' 
#' @return An object of class 'bvecmodel'.
#' 
#' @export
#' @method scale_error_correction bvecmodel
scale_error_correction.bvecmodel <- function(object, ...) {
  
  w <- object[["data"]][["train"]][["w"]]

  # Scaling twice is a no-op on the numbers -- dividing by the standard
  # deviation of the differences makes that standard deviation one, so the
  # second pass divides by one -- but it recomputes the factors from the series
  # it is given and overwrites the stored ones with those ones. The first
  # scaling would then be irreversible. Refused, the way a second rescaling is
  # by the attribute being dropped.
  if (!is.null(attr(w, "scale"))) {
    stop("The series in the error correction term are already scaled. Use ",
         "'rescale_error_correction' to put them back on the scale of the ",
         "input data first.")
  }

  # A cointegration space prior other than a multiple of the identity has a
  # direction, and that direction is stated in the units of 'w' it was built
  # against -- for coint$p_tau_i = "ml" the maximum likelihood estimate on the
  # series as they were. Scaling afterwards would leave it pointing somewhere
  # else without any sign of it. That holds for the constant model's p_tau_inv
  # and for the time varying model's transition p_tau alike.
  has_direction <- function(m) {
    !is.null(m) && !isTRUE(all.equal(m, diag(m[1, 1], nrow(m)), check.attributes = FALSE))
  }
  if (has_direction(object[["priors"]][["beta"]][["p_tau_inv"]]) ||
      has_direction(object[["priors"]][["beta"]][["p_tau"]])) {
    stop("The model already has a cointegration space prior that depends on the ",
         "scale of the error correction term. Call 'scale_error_correction' ",
         "before 'add_priors'.")
  }

  tt <- nrow(w)
  
  rescale_factors <- rep(1, ncol(w))
  pos_non_trend <- which(!dimnames(w)[[2]] %in% c("const", "trend"))
  rescale_factors[pos_non_trend] <- apply(diff(w[, pos_non_trend]), 2, sd)
  
  if ("trend" %in% dimnames(w)[[2]]) {
    pos_trend <- which(dimnames(w)[[2]] == "trend")
    rescale_factors[pos_trend] <- sd(w[, "trend"])
  }
  
  names(rescale_factors) <- dimnames(w)[[2]]
  attr(object[["data"]][["train"]][["w"]], "scale") <- rescale_factors
  
  w <- w / t(matrix(rescale_factors, length(rescale_factors), tt))
  
  object[["data"]][["train"]][["w"]][] <- w[]
  
  return(object)
}