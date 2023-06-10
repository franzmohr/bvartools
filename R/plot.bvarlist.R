#' Plotting Posterior Draws of Bayesian VAR or VEC Models
#' 
#' A plot function for objects of class \code{"bvarlist"}.
#' 
#' @param x an object of class \code{"bvarlist"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param model indicates, for which model plots should be produced. Can be a character
#' vector of country names or an integer vector of the positions of the elements in
#' argument \code{x}.
#' @param ... further graphical parameters.
#' 
#' @export 
plot.bvarlist <- function(x, ci = 0.95, type = "hist", model = NULL, ...) {
  
  if (is.null(model)) {
    for (i in 1:length(x)) {
      plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...)
    }
  } else {
    if (any(class(model) %in% c("numeric", "integer"))) {
      for (i in model) {
        plot(x[[i]], ci = ci, type = type, ctry = names(x)[i], ...)
      }
    }
  }
}


