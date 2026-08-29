#' Adding Error Bars to a Plot
#'
#' Adds the credible band of a statistic to an existing plot, which is marked by
#' the median and the mean of its posterior draws.
#'
#' @param pos a vector of positions of the error bars on the axis, along which they
#' are arranged. This is the x-axis for vertical and the y-axis for horizontal
#' error bars.
#' @param values a data frame or matrix with columns \code{"mean"}, \code{"median"},
#' \code{"qlower"} and \code{"qupper"} and one row per element of \code{pos}.
#' @param cap half of the width of the lines, which mark the bounds and the median
#' of an error bar.
#' @param labels a vector of labels, which are added at the upper end of the error
#' bars. If \code{NULL}, no labels are added.
#' @param col a vector of colours, which is recycled over the error bars.
#' @param lwd a vector of line widths, which is recycled over the error bars.
#' @param pch the plotting symbol used for the mean.
#' @param cex the size of the plotting symbol used for the mean.
#' @param cex_lab the size of the labels at the upper end of the error bars.
#' @param horizontal logical. If \code{TRUE}, the error bars are plotted
#' horizontally, i.e. the statistic is measured on the x-axis. Defaults to
#' \code{FALSE}.
#'
#' @noRd
draw_error_bars <- function(pos, values, cap, labels = NULL, col = "black", lwd = 1,
                            pch = 20, cex = 1, cex_lab = .8, horizontal = FALSE) {

  # Lower bound, median and upper bound. Colours and line widths are recycled
  # in the same order as the positions
  bounds <- unlist(values[, c("qlower", "median", "qupper")], use.names = FALSE)

  if (horizontal) {
    # Credible band
    graphics::segments(values[, "qlower"], pos, values[, "qupper"], pos,
                       col = col, lwd = lwd)
    graphics::segments(bounds, rep(pos - cap, 3), bounds, rep(pos + cap, 3),
                       col = col, lwd = lwd)
    # Mean
    graphics::points(values[, "mean"], pos, pch = pch, cex = cex, col = col)

    if (!is.null(labels)) {
      graphics::text(values[, "qupper"], pos, labels = labels, pos = 4, offset = .3,
                     cex = cex_lab, col = col)
    }
  } else {
    # Credible band
    graphics::segments(pos, values[, "qlower"], pos, values[, "qupper"],
                       col = col, lwd = lwd)
    graphics::segments(rep(pos - cap, 3), bounds, rep(pos + cap, 3), bounds,
                       col = col, lwd = lwd)
    # Mean
    graphics::points(pos, values[, "mean"], pch = pch, cex = cex, col = col)

    if (!is.null(labels)) {
      graphics::text(pos, values[, "qupper"], labels = labels, pos = 3, offset = .3,
                     cex = cex_lab, col = col)
    }
  }
}
