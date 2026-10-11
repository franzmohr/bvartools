#' Plotting Forecasts of a List of Bayesian VAR Models
#'
#' A plot function for objects of class 'expandwindbvarprdlist'.
#'
#' @param x an object of class 'expandwindbvarprdlist', usually, a result of a call to
#' \code{\link{predict.expandingwindow}}.
#' @param n_pre number of plotted observations that precede the forecasts.
#' If \code{NULL} (default), all available observations will be plotted.
#' @param ci interval used to calculate the credible bands of the forecasts.
#' @param ... further graphical parameters, which are passed on to
#' \code{\link[stats]{plot.ts}}.
#' 
#' @return \code{x}, invisibly. The function is called for its side effect, the
#' plot of the forecasts.
#'
#' @export 
#' @method plot expandwindbvarprdlist
plot.expandwindbvarprdlist <- function(x, n_pre = NULL, ci = 0.95, ...) {
  # The forecasts of the windows are a list of forecasts like any other, and
  # this used to be a line-for-line copy of bayests' method for one, which the
  # two would have drifted apart from.
  as_list <- x
  class(as_list) <- c("bvarprdlist", setdiff(class(x), "expandwindbvarprdlist"))
  plot(as_list, n_pre = n_pre, ci = ci, ...)
  invisible(x)
}
