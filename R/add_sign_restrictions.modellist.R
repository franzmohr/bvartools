#' @include add_sign_restrictions.R
NULL

#' @rdname add_sign_restrictions.bvarmodel
#' @export
#' @method add_sign_restrictions modellist
add_sign_restrictions.modellist <- function(object, ...) {

  object <- lapply(object, add_sign_restrictions, ...)
  class(object) <- list("modellist", "list")

  return(object)
}
