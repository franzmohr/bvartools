
#' @export
get_model_specifications.selcrit <- function(object, ...) {

  # The class attribute of a 'selcrit' object contains the classes of the model,
  # for which the selection criteria were calculated. Dispatch on the classes
  # that follow 'selcrit' finds the method of the model class, which may also be
  # defined in another package
  classes <- class(object)
  classes <- classes[-seq_len(match("selcrit", classes))]

  if (!any(classes != "list")) {
    stop("The class attribute of argument 'object' does not contain the class of its model.")
  }

  class(object) <- classes

  return(get_model_specifications(object, ...))
}
