#' @include historical_decomposition.R map_collection.R
NULL

#' @rdname historical_decomposition.bvarmodel
#' @export
#' @method historical_decomposition modellist
historical_decomposition.modellist <- function(x, ...) {
  .analyse_collection(x, historical_decomposition, "historical_decomposition", ...)
}
