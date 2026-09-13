#' Choose Best Model
#' 
#' A generic function that applies model selection criteria to find the best
#' model.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{choose_best_model.selcritlist}}.
#'
#' @export
choose_best_model <- function (object, ...) {
  UseMethod("choose_best_model")
}
