#' Transform a VEC Model to a VAR in Levels
#'
#' A generic function used to transform a vector error correction model into
#' its VAR form. The function invokes particular methods which depend on
#' the class of the first argument.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{vec_to_var.bvecmodel}},
#' \code{\link{vec_to_var.modellist}}.
#'
#' @export
vec_to_var <- function (object, ...) {
  UseMethod("vec_to_var")
}
