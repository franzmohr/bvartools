#' Summarising Bayesian Models
#'
#' summary method for class \code{"posteriorlist"}.
#'
#' @param object an object of class \code{"posteriorlist"}, usually, a result of a call to
#' \code{\link{draw_posterior}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list of class \code{"summary.posteriorlist"}.
#'
#' @export
summary.posteriorlist <- function(object, ...){
  
  object <- lapply(object, summary, ...)
  
  return(object)
}