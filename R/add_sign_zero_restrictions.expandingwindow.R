#' @include add_sign_zero_restrictions.R
NULL

#' @rdname add_sign_zero_restrictions.bvarmodel
#' @export
#' @method add_sign_zero_restrictions expandingwindow
add_sign_zero_restrictions.expandingwindow <- function(object, ...) {

  object <- lapply(object, add_sign_zero_restrictions, ...)
  class(object) <- list("expandingwindow", "list")

  return(object)
}
