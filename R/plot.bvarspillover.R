#' Plot a Spillover Index
#'
#' Plots the directional connectedness measures of Diebold and Yilmaz (2012).
#'
#' @param x an object of class 'bvarspillover', usually the result of a call to
#' \code{\link{spillover}}.
#' @param which which measure to plot. One of \code{"net"} (default),
#' \code{"to"} or \code{"from"}.
#' @param ... further arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @details One bar per variable, with the credible interval drawn on it where
#' the object carries one. Net values above zero mark variables that transmit
#' more forecast error variance than they receive.
#'
#' @return \code{NULL}, invisibly. Called for its side effect.
#'
#' @export
plot.bvarspillover <- function(x, which = "net", ...) {

  which <- match.arg(which, c("net", "to", "from"))
  values <- x[[which]]

  if (is.matrix(values) && !is.null(rownames(values))) {
    mid <- values[rownames(values) == "50%", ]
    low <- values[1, ]
    high <- values[nrow(values), ]
  } else {
    # keep_draws = TRUE: summarise here rather than refusing to plot.
    mid <- apply(values, 2, stats::median)
    ci <- x[["specification"]][["ci"]]
    low <- apply(values, 2, stats::quantile, probs = (1 - ci) / 2)
    high <- apply(values, 2, stats::quantile, probs = 1 - (1 - ci) / 2)
  }

  main <- switch(which,
                 "net" = "Net spillovers",
                 "to" = "Spillovers to others",
                 "from" = "Spillovers from others")

  dots <- list(...)
  args <- list(ylab = "Percentage", xlab = "", main = main,
               ylim = range(c(0, low, high)) * 1.1)
  args <- args[!(names(args) %in% names(dots))]

  at <- do.call(graphics::barplot, c(list(mid), args, dots))

  # The interval, on the bars rather than beside them.
  graphics::arrows(at, low, at, high, angle = 90, code = 3,
                   length = 0.05, col = "grey30")
  graphics::abline(h = 0)

  invisible(NULL)
}
