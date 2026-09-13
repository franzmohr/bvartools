#' Transform a VEC Model to a VAR in Levels
#'
#' The elements of an object of class \code{'modellist'} are transformed into
#' their VAR representation in levels.
#'
#' @param object an object of class \code{'modellist'}.
#' @param ... arguments passed forward to method.
#'
#' @return An object of class \code{'modellist'}.
#'
#' @export
#' @method vec_to_var modellist
vec_to_var.modellist <- function(object, ...) {

  class_of_object <- class(object)

  object <- lapply(object, vec_to_var, ...)

  class(object) <- class_of_object

  return(object)
}



#' Transform VEC Models to VARs in Levels
#'
#' The models of the estimation windows of an object of class
#' \code{'expandingwindow'} are transformed into their VAR representation in
#' levels, which is the form in which forecasts, forecast errors and
#' out-of-sample selection criteria are obtained for VEC models.
#'
#' @param object an object of class \code{'expandingwindow'}.
#' @param ... arguments passed forward to method.
#'
#' @return An object of class \code{'expandingwindow'}.
#'
#' @export
#' @method vec_to_var expandingwindow
vec_to_var.expandingwindow <- function(object, ...) {

  # Without this method an expanding window of VEC models had no way to its
  # forecasts: add_forecast_input() refuses a 'bvecmodel' and points to
  # vec_to_var(), which in turn had no method for the windows.
  class_of_object <- class(object)

  object <- lapply(object, vec_to_var, ...)

  class(object) <- class_of_object

  return(object)
}