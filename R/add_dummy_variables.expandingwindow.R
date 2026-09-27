#' @include add_dummy_variables.R
NULL

#' @rdname add_dummy_variables.bvarmodel
#' @export
#' @method add_dummy_variables expandingwindow
add_dummy_variables.expandingwindow <- function(object, impulse = NULL, step = NULL, data = NULL, ...) {

  # A window that ends before the period of a dummy leaves it out rather than
  # refusing it, since a forecast made at its end could not have known about
  # the event. The last window has the whole sample, so a dummy it refuses
  # would have been refused by every window.
  n_windows <- length(object)
  for (i in seq_len(n_windows)) {
    frequency <- stats::frequency(object[[i]][["data"]][["train"]][["y"]])
    dummies <- .dummy_specification(impulse, step, data, frequency)
    empty <- if (i == n_windows) "stop" else "drop"
    object[[i]] <- .add_dummy_variables(object[[i]], dummies, empty = empty)
  }
  class(object) <- list("expandingwindow", "list")

  return(object)
}
