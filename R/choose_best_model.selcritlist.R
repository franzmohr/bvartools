#' Plotting Selection Criteria
#' 
#' A plot function for objects of class \code{"selcritlist"}.
#' 
#' @param x object object of class \code{"selcritlist"}, usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default), \code{"HQ"}.
#' 
#' @export
choose_best_model.selcritlist <- function(object, criterion = "BIC", ...) {
  
  res <- lapply(object, function(y) {y[["summary"]]})
  res <- do.call("rbind", res)
  
  pos <- which(res[, criterion] == max(res[, criterion]))
  
  return(pos)
}
