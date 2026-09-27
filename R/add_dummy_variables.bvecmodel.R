#' @include add_dummy_variables.R
NULL

#' @rdname add_dummy_variables.bvarmodel
#' @export
#' @method add_dummy_variables bvecmodel
add_dummy_variables.bvecmodel <- function(object, impulse = NULL, step = NULL, data = NULL, ...) {

  frequency <- stats::frequency(object[["data"]][["train"]][["y"]])
  dummies <- .dummy_specification(impulse, step, data, frequency)

  return(.add_dummy_variables(object, dummies))
}
