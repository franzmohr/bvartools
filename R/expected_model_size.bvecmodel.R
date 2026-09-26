#' @include expected_model_size.R
#' @rdname expected_model_size.bvarmodel
#' @export
expected_model_size.bvecmodel <- function(object, chains = NULL, ...) {
  .model_size(object, chains = chains)
}
