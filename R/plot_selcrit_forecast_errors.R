#' Plotting Forecast Errors of Selection Criteria Objects
#'
#' Plots the out-of-sample statistics of an object of class 'selcritlist' as
#' error bars for each forecast horizon and endogenous variable.
#'
#' @param x an object of class 'selcritlist'.
#' @param criterion the out-of-sample statistic that should be plotted. Available
#' choices are \code{"FE"}, \code{"AFE"} and \code{"RSFE"}.
#' @param col a vector of colours, which is recycled over the models in \code{x}.
#' @param pch the plotting symbol used for the mean of a model.
#' @param cex the size of the plotting symbol used for the mean of a model.
#' @param lwd a vector of line widths, which is recycled over the models in \code{x}.
#' @param main the title of the plot. If \code{NULL} (default), no title is added.
#' @param xlab the label of the x-axis.
#' @param ... further graphical parameters, which are passed on to
#' \code{\link[graphics]{plot}}.
#'
#' @details Each row of the plot corresponds to an endogenous variable and the
#' x-axis to the forecast horizon. Within a forecast horizon each model is
#' represented by an error bar, which covers the credible band of the respective
#' statistic and which is marked by the median and the mean of its posterior draws.
#' The position of a model in \code{x} is added above its error bar. For
#' \code{criterion = "FE"} a horizontal reference line is added at zero.
#'
#' @noRd
plot_selcrit_forecast_errors <- function(x, criterion, col = "black", pch = 20,
                                        cex = 1, lwd = 1, main = NULL, xlab = "Forecast horizon", ...) {

  n_models <- length(x)
  col <- rep_len(col, n_models)
  lwd <- rep_len(lwd, n_models)

  errors <- lapply(x, function(y, criterion) {y[[criterion]]}, criterion = criterion)
  # Models, which do not contain the statistic, are omitted from the calculations
  # below, but their positions in 'x' are maintained in object 'errors'
  avail <- errors[!unlist(lapply(errors, is.null))]

  vars <- unique(unlist(lapply(avail, function(y) {as.character(y[, "variable"])})))
  n_vars <- length(vars)
  n_ahead <- sort(unique(unlist(lapply(avail, function(y) {y[, "h"]}))))

  # Horizontal positions of the error bars of the models within a forecast horizon.
  # The width must be clearly smaller than one so that the error bars of a
  # forecast horizon can be distinguished from those of the neighbouring horizons.
  bar_width <- .6
  if (n_models == 1) {
    offset <- 0
  } else {
    offset <- seq(from = -bar_width / 2, to = bar_width / 2, length.out = n_models)
  }
  # Forecast errors can be positive and negative, so a reference line is added
  zero_line <- criterion == "FE"

  # Common limits of the axes
  xlim <- range(n_ahead) + c(-.5, .5)

  # Half of the width of the horizontal lines of an error bar
  cap <- ifelse(n_models == 1, .15, min(.15, bar_width / (n_models - 1) / 3))
  cap <- min(cap, diff(xlim) / 20)

  ylim <- lapply(vars, function(v, avail, zero_line) {
    temp <- unlist(lapply(avail, function(y, v) {
      as.matrix(y[y[, "variable"] == v, c("mean", "median", "qlower", "qupper")])
    }, v = v))
    temp <- temp[is.finite(temp)]
    if (zero_line) {
      # The reference line must be visible
      temp <- c(temp, 0)
    }
    temp <- range(temp)
    if (diff(temp) == 0) {
      temp <- temp + c(-.5, .5)
    }
    # Space for the model numbers above the error bars
    temp + c(0, .1 * diff(temp))
  }, avail = avail, zero_line = zero_line)
  names(ylim) <- vars

  # 'layout' is called below, so all parameters have to be restored on exit
  orig_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(orig_par))

  mar_ax_x <- 3.5
  mar_ax_y <- 3
  mar_min <- .5
  # Space above the subplots is only required, if a title is added
  oma_top <- ifelse(is.null(main), 0, 2.5)

  graphics::par(oma = c(0, 0, oma_top, 0))

  # Space that is required for the axes, measured in fractions of the region
  # that is available for the subplots
  csi <- graphics::par("csi")
  ax_height <- mar_ax_x * csi / (graphics::par("din")[2] - oma_top * csi)

  lab_size <- .05
  plot_height <- (1 - ax_height) / n_vars

  layout_matrix <- cbind(1:n_vars, n_vars + 1:n_vars)
  graphics::layout(layout_matrix,
                   widths = c(lab_size, 1 - lab_size),
                   heights = c(rep(plot_height, n_vars - 1), ax_height + plot_height))

  # Variable names left of the rows of the plot
  for (i in 1:n_vars) {
    graphics::par(mar = c(ifelse(i == n_vars, mar_ax_x, mar_min), 0, mar_min, 0))
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = vars[i], adj = 0.5, srt = 90)
  }

  # The x-axis is not extended beyond 'xlim' so that the distance between the
  # outer error bars and the borders of a subplot is the same as the distance
  # between the error bars and the separators of the forecast horizons
  graphics::par(xaxs = "i")

  for (i in 1:n_vars) {
    graphics::par(mar = c(ifelse(i == n_vars, mar_ax_x, mar_min),
                          mar_ax_y, mar_min, mar_min))

    graphics::plot(NA, xlim = xlim, ylim = ylim[[vars[i]]], axes = FALSE,
                   main = "", xlab = "", ylab = "", ...)

    # Separators between the forecast horizons
    graphics::abline(v = (n_ahead[-length(n_ahead)] + n_ahead[-1]) / 2, lty = 3, col = "darkgrey")

    if (zero_line) {
      graphics::abline(h = 0, col = "darkgrey")
    }

    for (j in 1:n_models) {
      temp <- errors[[j]]
      if (is.null(temp)) {
        # Statistic is not contained in the respective model
        next
      }
      temp <- temp[temp[, "variable"] == vars[i], ]
      if (nrow(temp) == 0) {
        # Variable is not contained in the respective model
        next
      }

      pos <- temp[, "h"] + offset[j]

      # The label of an error bar is the position of the model in 'x'
      draw_error_bars(pos, temp, cap = cap, labels = j,
                      col = col[j], lwd = lwd[j], pch = pch, cex = cex)
    }

    if (i == n_vars) {
      graphics::axis(1, at = n_ahead)
      if (!is.null(xlab)) {
        graphics::mtext(xlab, side = 1, line = 2.2)
      }
    }
    graphics::axis(2)
    graphics::box()
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = .5,
                    cex = graphics::par("cex.main") / graphics::par("cex"),
                    font = graphics::par("font.main"))
  }
}
