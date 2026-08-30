#' Summarising Bayesian Models
#'
#' summary method for class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' The method passes the last element of the expanding window model list to its
#' corresponding summary function; assuming that this object contains the model
#' with the highest amount of available observations.
#' 
#'
#' @return A list of class 'summary.bvarmodel' or 'summary.bvecmodel'.
#'
#' @export
summary.expandingwindow <- function(object, ...){
  
  # Only evaluate the last model with the highest amount of observations
  object <- summary(object[[length(object)]])
  
  return(object)
}