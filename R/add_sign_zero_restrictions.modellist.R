#' @include add_sign_zero_restrictions.R
NULL

#' @rdname add_sign_zero_restrictions.bvarmodel
#' @export
#' @method add_sign_zero_restrictions modellist
add_sign_zero_restrictions.modellist <- function(object, ...) {

  object <- lapply(object, add_sign_zero_restrictions, ...)
  class(object) <- list("modellist", "list")

  return(object)
}
