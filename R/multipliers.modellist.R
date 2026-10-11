#' @include multipliers.R map_collection.R
NULL

#' @rdname multipliers.bvarmodel
#' @export
#' @method multipliers modellist
multipliers.modellist <- function(x, ...) {
  .analyse_collection(x, multipliers, "multipliers", ...)
}
