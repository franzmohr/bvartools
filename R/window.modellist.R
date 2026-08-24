#' Time Series Windows
#' 
#' Forwards model input to the same function for individual models.
#' 
#' @param x an object of class 'modellist'.
#' @param start the start time of the period of interest.
#' @param end the end time of the period of interest.
#' 
#' @return An object of class 'modellist'.
#' 
#' @export
#' @method window modellist
window.modellist <- function(x, start = NULL, end = NULL, ...) {
  
  for (i in 1:length(x)) {
    x[[i]] <- stats::window(x = x[[i]], start = start, end = end)
  }
  
  return(x)
}