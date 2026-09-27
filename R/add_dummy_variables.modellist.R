#' @include add_dummy_variables.R
NULL

#' @rdname add_dummy_variables.bvarmodel
#' @export
#' @method add_dummy_variables modellist
add_dummy_variables.modellist <- function(object, ...) {

  object <- lapply(object, add_dummy_variables, ...)
  class(object) <- list("modellist", "list")

  return(object)
}
