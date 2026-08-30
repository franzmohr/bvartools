#' Summarising Bayesian Models
#'
#' summary method for class 'modellist'.
#'
#' @param object an object of class 'modellist'.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list of class 'summary.modellist'.
#'
#' @export
#' @method summary modellist
summary.modellist <- function(object, ...){
  
  object <- lapply(object, summary, ...)
  
  class(object) <- c("summary.modellist", "list")
  
  return(object)
}