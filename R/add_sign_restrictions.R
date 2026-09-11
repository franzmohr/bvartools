#' Sign Restrictions
#' 
#' A generic function used to identify the shocks of a model by the signs of
#' the impulse responses they produce.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
add_sign_restrictions <- function(object, ...) {
  UseMethod("add_sign_restrictions")
}
