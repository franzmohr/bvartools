#' @include ssvs_prior.R map_collection.R
NULL

#' @rdname ssvs_prior.bvarmodel
#' @export
#' @method ssvs_prior modellist
ssvs_prior.modellist <- function(object, ...) {
  .map_collection(object, ssvs_prior, ...)
}
