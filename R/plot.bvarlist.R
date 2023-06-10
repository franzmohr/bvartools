#' Plotting Posterior Draws of Bayesian VAR or VEC Models
#' 
#' A plot function for objects of class \code{"bvarlist"}.
#' 
#' @param x an object of class \code{"bvarlist"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param model numeric or integer indicating for which models in argument \code{"x"} plots should be produced.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.bvarlist <- function(x, ci = 0.95, type = "hist", model = NULL, ...) {
  
  if (is.null(model)) {
    model <- 1:length(x)
  } else {
    if (!any(c("numeric", "integer") %in% class(model))) {
      stop("If specified, argument 'model' must be numeric or integer.")
    }
  }
  
  for (i in model) {
    if (!is.null(x[[i]][["error"]])) {
      if (x[[i]][["error"]]) {
        next
      }
    } else {
      plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...) 
    }
  }
}


