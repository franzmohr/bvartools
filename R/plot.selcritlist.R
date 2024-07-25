#' Plotting Selection Criteria
#' 
#' A plot function for objects of class \code{"selcritlist"}.
#' 
#' @param x an object of class \code{"selcritlist"}, usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"}, \code{"HQ"}.
#' @param ... further graphical parameters passed on to \link[graphics]{boxplot}.
#' 
#' @export
plot.selcritlist <- function(x, criterion = "LL", ...) {
  
  if (criterion == "LL") {
    res <- lapply(x, function(x) {rowSums(x$ll)})
  }
  
  if (criterion == "AIC") {
    res <- lapply(x, function(x) {2 * x$npara - 2 * rowSums(x$ll)})
  }
  
  if (criterion == "BIC") {
    res <- lapply(x, function(x) {log(ncol(x$ll)) * x$npara - 2 * rowSums(x$ll)})
  }
  
  if (criterion == "HQ") {
    res <- lapply(x, function(x) {log(log(ncol(x$ll))) * x$npara - 2 * rowSums(x$ll)})
  }
  
  # Model specifications for x-axis
  model_specs <- NULL
  for (i in 1:length(x)) {
    names_i <- names(x[[i]][["summary"]])
    pos_i <- 1:(which(names_i == "LL") - 1)
    names_i <- names_i[pos_i]
    model_specs <- append(model_specs, paste0(paste0(names_i, " = ", x[[i]][["summary"]][,pos_i]), collapse = "\n"))
  }
  names(res) <- model_specs

  graphics::boxplot(res, ...)
}
