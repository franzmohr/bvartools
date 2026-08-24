#' Time Series Windows
#' 
#' Restricts the observations in the training sample of a model of class
#' 'bvarmodel' to the specified time window.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param start the start time of the period of interest.
#' @param end the end time of the period of interest.
#' 
#' @return An object of class 'bvarmodel'.
#' 
#' @export
#' @method window bvarmodel
window.bvarmodel <- function(x, start = NULL, end = NULL, ...) {
  
  
  if (!is.null(x[["data"]][["train"]][["z"]])) {
    k <- x[["model"]][["k"]]
    orig_time <- stats::time(x[["data"]][["train"]][["y"]])
  }
  
  x[["data"]][["train"]][["y"]] <- stats::window(x[["data"]][["train"]][["y"]],
                                                         start = start, end = end, ...)
  
  if (!is.null(x[["data"]][["train"]][["x"]])) {
    x[["data"]][["train"]][["x"]] <- stats::window(x[["data"]][["train"]][["x"]],
                                                        start = start, end = end, ...) 
  }
  
  if (!is.null(x[["data"]][["train"]][["z"]])) {
    new_time <- stats::time(x[["data"]][["train"]][["y"]])
    pos <- which(orig_time %in% new_time)
    pos <- rep(pos * k, each = k) - k + 1 + rep(0:(k - 1), length(pos))
    x[["data"]][["train"]][["z"]] <- x[["data"]][["train"]][["z"]][pos,]
  }
  
  return(x)
}