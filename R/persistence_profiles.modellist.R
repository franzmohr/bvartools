#' @include persistence_profiles.R
NULL

#' @rdname persistence_profiles.bvecmodel
#' @export
#' @method persistence_profiles modellist
persistence_profiles.modellist <- function(object, ...) {
  object <- lapply(object, persistence_profiles, ...)
  class(object) <- list("modellist", "list")
  object
}

#' @rdname persistence_profiles.bvecmodel
#' @export
#' @method persistence_profiles expandingwindow
persistence_profiles.expandingwindow <- function(object, ...) {
  object <- lapply(object, persistence_profiles, ...)
  class(object) <- list("expandingwindow", "list")
  object
}
