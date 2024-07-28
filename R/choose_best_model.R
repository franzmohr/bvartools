#' Choose Best Model
#' 
#' A generic function that applies model selection criteria to find the best
#' model.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
choose_best_model <- function (object, ...) {
  UseMethod("choose_best_model")
}