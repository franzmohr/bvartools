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