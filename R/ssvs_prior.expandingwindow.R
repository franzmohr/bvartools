#' @include ssvs_prior.R map_collection.R
NULL

#' @rdname ssvs_prior.bvarmodel
#' @export
#' @method ssvs_prior expandingwindow
ssvs_prior.expandingwindow <- function(object, ...) {
  .map_collection(object, ssvs_prior, ...)
}
