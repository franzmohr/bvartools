#' @include expected_size.R
#' @rdname expected_size.bvarmodel
#' @export
expected_size.bvecmodel <- function(object, chains = NULL, ...) {
  .model_size(object, chains = chains)
}
