#' @include expected_size.R
#' @rdname expected_size.bvarmodel
#' @export
expected_size.modellist <- function(object, ...) {
  .model_list_size(object, ...)
}
