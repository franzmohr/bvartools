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