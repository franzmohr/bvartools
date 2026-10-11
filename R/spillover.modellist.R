#' Spillover Index
#'
#' Produces the connectedness measures of Diebold and Yilmaz (2012) for every
#' model in an object of class 'modellist'.
#'
#' @param object an object of class 'modellist'.
#' @param ... arguments passed forward to \code{\link{spillover.bvarmodel}}.
#'
#' @return A list of objects of class 'bvarspillover', one per model. An error
#' in one model stops the call and names the model.
#'
#' @export
spillover.modellist <- function(object, ...) {
  # Every error used to become a NULL element, under a guard for an 'error'
  # flag nothing sets: a misspelt argument returned a list of NULLs without a
  # word.
  .analyse_collection(object, spillover, "spillover", ...)
}
