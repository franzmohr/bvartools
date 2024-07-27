#' Choose Best Model
#' 
#' Chooses the best model according the selection criteria in an object of class
#' 'selcritlist'.
#' 
#' @param object object of class 'selcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default), \code{"HQ"}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' If argument \code{criterion} is "LL", the model with the maximum value is chosen,
#' otherwise, the model with the minimum value.
#' 
#' @returns An integer giving the position of the best model in the list provided
#' in argument \code{object}.
#' 
#' @export
choose_best_model.selcritlist <- function(object, criterion = "BIC", ...) {
  
  res <- lapply(object, function(y) {y[["summary"]]})
  res <- do.call("rbind", res)
  
  if (criterion == "LL") {
    pos <- which(res[, criterion] == max(res[, criterion]))
  } else {
    pos <- which(res[, criterion] == min(res[, criterion])) 
  }
  
  return(pos)
}
