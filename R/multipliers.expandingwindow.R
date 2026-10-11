#' @include multipliers.R map_collection.R
NULL

#' @rdname multipliers.bvarmodel
#' @export
#' @method multipliers expandingwindow
multipliers.expandingwindow <- function(x, ...) {
  .analyse_collection(x, multipliers, "multipliers", ...)
}
