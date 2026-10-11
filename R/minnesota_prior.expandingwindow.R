#' @include minnesota_prior.R map_collection.R
NULL

#' @rdname minnesota_prior.modellist
#' @export
#' @method minnesota_prior expandingwindow
minnesota_prior.expandingwindow <- function(object, ...) {
  .map_collection(object, minnesota_prior, ...)
}
