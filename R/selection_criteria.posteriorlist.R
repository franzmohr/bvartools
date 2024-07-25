#' Model Selection Criteria
#'
#' Calculates model selection criteria for a list of Bayesian models.
#'
#' @param object an object of class \code{"posteriorlist"}, usually, the
#' output of a call to \code{\link{draw_posterior}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list.
#' 
#' @export
selection_criteria.posteriorlist <- function(object, ...){
  
  object <- lapply(object, selection_criteria)
  
  class(object) <- append("selcritlist", class(object))
  
  return(object)
}