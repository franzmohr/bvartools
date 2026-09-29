#' @include time_variation_test.R
NULL

#' @rdname time_variation_test.bvarmodel
#' @export
#' @method time_variation_test modellist
time_variation_test.modellist <- function(object, ...) {
  object <- lapply(object, time_variation_test, ...)
  class(object) <- list("modellist", "list")
  object
}

#' @rdname time_variation_test.bvarmodel
#' @export
#' @method time_variation_test expandingwindow
time_variation_test.expandingwindow <- function(object, ...) {
  object <- lapply(object, time_variation_test, ...)
  class(object) <- list("expandingwindow", "list")
  object
}
