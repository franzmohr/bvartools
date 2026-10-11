#' @include inclusion_prior.R map_collection.R
NULL

#' @rdname inclusion_prior.bvarmodel
#' @export
#' @method inclusion_prior modellist
inclusion_prior.modellist <- function(object, ...) {
  .map_collection(object, inclusion_prior, ...)
}
