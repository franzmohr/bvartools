#' Model Selection Criteria
#'  
#' A generic function used to calculate selection criteria.
#' 
#' @param object an object of class \code{"posteriorlist"}, \code{"bvarmodel"} or
#' \code{"bvecmodel"}.
#' @param ... arguments passed forward to method.
#' 
#' @export
selection_criteria <- function (object, ...) {
 UseMethod("selection_criteria")
}