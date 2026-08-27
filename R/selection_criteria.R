#' Model Selection Criteria
#' 
#' Generic function used to calculate selection criteria.
#' 
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @export
selection_criteria <- function (object, ...) {
 UseMethod("selection_criteria")
}