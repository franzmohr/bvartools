#' Plotting Forecasts of a List of Bayesian VAR Models
#'
#' A plot function for objects of class 'bvarprdlist'.
#'
#' @param x an object of class 'bvarprdlist', usually, a result of a call to
#' \code{\link{predict.modellist}}.
#' @param n_pre number of plotted observations that precede the forecasts.
#' If \code{NULL} (default), all available observations will be plotted.
#' @param ci interval used to calculate the credible bands of the forecasts.
#' @param ... further graphical parameters, which are passed on to
#' \code{\link[stats]{plot.ts}}.
#'
#' @details The forecasts of the models in \code{x} are plotted on a grid, where
#' each column corresponds to a model and each row to an endogenous variable.
#'
#' @export
#' @method plot bvarprdlist
plot.bvarprdlist <- function(x, n_pre = NULL, ci = 0.95, ...) {

  if (!all(unlist(lapply(x, function(x) {"bvarprd" %in% class(x)})))) {
    for (i in 1:length(x)) {
      plot(x[[i]], n_pre = n_pre, ci = ci, ...)
    }
    return(invisible(x))
  }

  dots <- list(...)

  # 'layout' is called below, so all parameters have to be restored on exit
  orig_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(orig_par))

  series <- lapply(x, forecast_plot_series, n_pre = n_pre, ci = ci)

  n_models <- length(series)
  model_names <- names(series)
  if (is.null(model_names)) {
    model_names <- as.character(1:n_models)
  }
  vars <- unique(unlist(lapply(series, names)))
  n_vars <- length(vars)

  # Common limits so that the axes of the subplots are aligned
  xlim <- range(unlist(lapply(series, function(x) {lapply(x, stats::time)})))
  ylim <- lapply(vars, function(i, series) {
    temp <- unlist(lapply(series, function(x, i) {x[[i]]}, i = i))
    range(temp[is.finite(temp)])
  }, series = series)
  names(ylim) <- vars

  # Space that is required for the axes of the subplots, measured in
  # fractions of the device size
  mar_ax <- 2.5
  mar_min <- .5
  ax_width <- mar_ax * graphics::par("csi") / graphics::par("din")[1]
  ax_height <- mar_ax * graphics::par("csi") / graphics::par("din")[2]

  lab_size <- .05
  plot_width <- (1 - lab_size - ax_width) / n_models
  plot_height <- (1 - lab_size - ax_height) / n_vars

  layout_matrix <- matrix(0, nrow = n_vars + 1, ncol = n_models + 1)
  layout_matrix[1, -1] <- 1:n_models
  layout_matrix[-1, 1] <- n_models + 1:n_vars
  layout_matrix[-1, -1] <- layout_matrix[n_vars + 1, 1] + 1:(n_models * n_vars)
  graphics::layout(layout_matrix,
                   widths = c(lab_size, ax_width + plot_width, rep(plot_width, n_models - 1)),
                   heights = c(lab_size, rep(plot_height, n_vars - 1), ax_height + plot_height))

  # Model names above the columns of the grid
  for (j in 1:n_models) {
    graphics::par(mar = c(0, ifelse(j == 1, mar_ax, mar_min), 0, mar_min))
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = model_names[j], adj = 0.5)
  }

  # Variable names left of the rows of the grid
  for (i in 1:n_vars) {
    graphics::par(mar = c(ifelse(i == n_vars, mar_ax, mar_min), 0, mar_min, 0))
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = vars[i], adj = 0.5, srt = 90)
  }

  # Subplots, drawn column-wise as specified in 'layout_matrix'
  for (j in 1:n_models) {
    for (i in 1:n_vars) {
      graphics::par(mar = c(ifelse(i == n_vars, mar_ax, mar_min),
                            ifelse(j == 1, mar_ax, mar_min),
                            mar_min, mar_min))

      if (is.null(series[[j]][[vars[i]]])) {
        # Variable is not contained in the respective model
        graphics::plot.new()
        next
      }

      # Defaults that may be overridden via '...'
      args <- list(plot.type = "single", lty = c(1, 2, 3, 2), main = "",
                   xlab = "", ylab = "", axes = FALSE,
                   xlim = xlim, ylim = ylim[[vars[i]]])
      args <- args[!(names(args) %in% names(dots))]

      do.call(stats::plot.ts, c(list(series[[j]][[vars[i]]]), args, dots))

      if (i == n_vars) {
        graphics::axis(1)
      }
      if (j == 1) {
        graphics::axis(2)
      }
      graphics::box()
    }
  }

  return(invisible(x))
}
