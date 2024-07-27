#' Choose Best Model
#' 
#' A generic function that applies model selection criteria to find the best
#' model.
#' 
#' @param object a list of estimated models.
#' @param ... arguments passed forward to method.
#' 
#' @export
choose_best_model <- function (object, ...) {
  UseMethod("choose_best_model")
}