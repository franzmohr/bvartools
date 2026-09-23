#' Sign and Zero Restrictions
#'
#' A generic function used to identify the shocks of a model by the signs of
#' the impulse responses they produce and by responses that are restricted to
#' be zero.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_sign_zero_restrictions.bvarmodel}}.
#'
#' @export
add_sign_zero_restrictions <- function(object, ...) {
  UseMethod("add_sign_zero_restrictions")
}
