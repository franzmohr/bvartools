#' Plotting Forecast Errors per Period
#'
#' Plots the forecast errors of a list of Bayesian models as error bars for each
#' period, for which a forecast was made.
#'
#' @param x an object of class 'bvarmodel', 'bvecmodel', 'expandingwindow' or
#' 'modellist', which contains forecast errors. Usually, a result of a call to
#' \code{\link{add_forecast_errors}}.
#' @param criterion the statistic that should be plotted. Available choices are
#' \code{"FE"} for forecast errors, \code{"AFE"} (default) for absolute forecast
#' errors and \code{"RSFE"} for root squared forecast errors.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible
#' band. Defaults to 0.95.
#' @param col a vector of colours, which is recycled over the models in \code{x}.
#' @param pch the plotting symbol used for the mean of a model.
#' @param cex the size of the plotting symbol used for the mean of a model.
#' @param lwd a vector of line widths, which is recycled over the models in \code{x}.
#' @param main the title of the plot. If \code{NULL} (default), no title is added.
#' @param xlab the label of the x-axis. If \code{NULL}, no label is added.
#' @param ... further graphical parameters, which are passed on to
#' \code{\link[graphics]{plot}}.
#'
#' @details Each row of the plot corresponds to an endogenous variable, each column
#' to a forecast horizon and the x-axis to the period, for which a forecast was
#' made. Within a period each model is represented by an error bar, which covers
#' the credible band of the respective statistic and which is marked by the median
#' and the mean of the posterior draws of its forecast errors. The position of a
#' model in \code{x} is added above its error bar. For \code{criterion = "FE"} a
#' horizontal reference line is added at zero.
#'
#' In contrast to \code{\link{plot.selcritlist}}, which plots statistics that are
#' calculated across all periods of the evaluation sample, this function plots the
#' statistics of each period separately. Thus, it is most useful for forecast errors,
#' which were obtained for multiple periods as, for example, in the case of
#' expanding window estimation. For objects of class 'expandingwindow' the forecast
#' errors of all estimation windows of a model are combined, where the period of an
#' error is obtained from the end of the training sample of the respective window.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' # AR(1) models as benchmark
#' model <- create_bvarmodel(data = us_macrodata,
#'                           p = 1:2,
#'                           deterministic = "none",
#'                           error = "gamma",
#'                           iterations = 10,
#'                           burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#'
#' # Obtain objects for expanding window estimation
#' model <- use_expanding_window(model, start = 2007)
#'
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(shape = 3, rate = 1))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#'
#' # Add forecasts
#' model <- add_posterior_forecasts(model)
#'
#' # Add forecast errors
#' model <- add_forecast_errors(model, test_sample = us_macrodata)
#'
#' # Plot absolute forecast errors of each period
#' plot_forecast_errors_by_period(model)
#'
#' @export
plot_forecast_errors_by_period <- function(x, criterion = "AFE", ci = .95,
                                           col = "black", pch = 20, cex = 1, lwd = 1,
                                           main = NULL, xlab = "Period", ...) {

  if (!criterion %in% c("FE", "AFE", "RSFE")) {
    stop("Argument 'criterion' must be either 'FE', 'AFE' or 'RSFE'.")
  }

  if (ci <= 0 | ci >= 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low

  if (!any(c("bvarmodel", "bvecmodel", "expandingwindow", "modellist") %in% class(x))) {
    stop("Argument 'x' must be an object of class 'bvarmodel', 'bvecmodel', 'expandingwindow' or 'modellist'.")
  }

  # Objects, which are not a list of models, contain a single model
  if (!"modellist" %in% class(x)) {
    x <- list(x)
  }

  n_models <- length(x)
  col <- rep_len(col, n_models)
  lwd <- rep_len(lwd, n_models)

  # Statistics of the posterior draws of the forecast errors of each variable,
  # period and forecast horizon
  errors <- lapply(x, function(y, criterion, probs) {

    # An object of class 'expandingwindow' contains one model per estimation window
    if (!"expandingwindow" %in% class(y)) {
      y <- list(y)
    }

    result <- NULL
    for (i in 1:length(y)) {

      values <- y[[i]][["posterior"]][["forecast_errors"]]
      if (is.null(values)) {
        # No test data were available for the respective window
        next
      }
      values <- as.matrix(values)

      varnames <- y[[i]][["model"]][["endogen"]]
      if (is.null(varnames)) {
        varnames <- dimnames(y[[i]][["data"]][["train"]][["y"]])[[2]]
      }
      k <- length(varnames)
      # The columns of the forecast errors are ordered by the endogenous variables
      # within the forecast horizons
      h <- rep(1:(ncol(values) / k), each = k)

      # Forecasts were made for the periods that follow the training sample
      tsp_train <- stats::tsp(y[[i]][["data"]][["train"]][["y"]])

      if (criterion == "AFE") {
        values <- abs(values)
      }
      if (criterion == "RSFE") {
        values <- values^2
      }
      temp <- t(apply(values, 2, function(z, probs) {
        z <- z[!is.na(z)]
        if (length(z) == 0) {
          return(rep(NA_real_, 3 + length(probs)))
        }
        c(mean(z), stats::median(z), stats::quantile(z, probs = probs))
      }, probs = probs))
      if (criterion == "RSFE") {
        temp <- sqrt(temp)
      }

      result <- rbind(result,
                      data.frame(variable = rep(varnames, max(h)),
                                 period = tsp_train[2] + h / tsp_train[3],
                                 h = h,
                                 mean = temp[, 1],
                                 median = temp[, 2],
                                 qlower = temp[, 3],
                                 qupper = temp[, 4],
                                 stringsAsFactors = FALSE))
    }

    if (!is.null(result)) {
      result <- result[!is.na(result[, "mean"]), ]
      if (nrow(result) == 0) {
        result <- NULL
      }
    }

    return(result)
  }, criterion = criterion, probs = c(ci_low, ci_high))

  # Models, which do not contain forecast errors, are omitted from the calculations
  # below, but their positions in 'x' are maintained in object 'errors'
  avail <- errors[!unlist(lapply(errors, is.null))]

  if (length(avail) == 0) {
    stop("Argument 'x' does not contain any forecast errors. You might want to use\nfunction 'add_forecast_errors' before this function.")
  }

  vars <- unique(unlist(lapply(avail, function(y) {y[, "variable"]})))
  n_vars <- length(vars)
  n_ahead <- sort(unique(unlist(lapply(avail, function(y) {y[, "h"]}))))
  # Periods are plotted at equidistant positions
  periods <- sort(unique(unlist(lapply(avail, function(y) {y[, "period"]}))))
  period_labels <- format(periods, trim = TRUE)

  # Horizontal positions of the error bars of the models within a period.
  # The width must be clearly smaller than one so that the error bars of a
  # period can be distinguished from those of the neighbouring periods.
  bar_width <- .6
  if (n_models == 1) {
    offset <- 0
  } else {
    offset <- seq(from = -bar_width / 2, to = bar_width / 2, length.out = n_models)
  }
  # Forecast errors can be positive and negative, so a reference line is added
  zero_line <- criterion == "FE"

  # Common limits of the axes
  xlim <- c(0, length(periods)) + .5

  # Half of the width of the horizontal lines of an error bar
  cap <- ifelse(n_models == 1, .15, min(.15, bar_width / (n_models - 1) / 3))
  cap <- min(cap, diff(xlim) / 20)

  # The limits of the y-axis are common across the forecast horizons of a variable
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

  mar_ax_y <- 3
  mar_min <- .5
  lab_size <- .05
  # Space above and below the subplots is only required, if a title or a label
  # of the x-axis is added
  oma_top <- ifelse(is.null(main), 0, 2.5)
  oma_bottom <- ifelse(is.null(xlab), 0, 1.8)

  graphics::par(oma = c(oma_bottom, 0, oma_top, 0))

  csi <- graphics::par("csi")
  in_width <- graphics::par("din")[1]
  in_height <- graphics::par("din")[2] - (oma_top + oma_bottom) * csi

  # The labels of the x-axis are turned by 90 degrees, so they extend downwards
  # by their width. It is measured before 'layout' reduces the size of the
  # plotted text and scaled with the respective factor below
  label_width <- max(graphics::strwidth(period_labels, units = "inches",
                                        cex = graphics::par("cex.axis")))

  layout_matrix <- matrix(0, nrow = n_vars + 1, ncol = length(n_ahead) + 1)
  layout_matrix[1, -1] <- 1:length(n_ahead)
  layout_matrix[-1, 1] <- length(n_ahead) + 1:n_vars
  layout_matrix[-1, -1] <- layout_matrix[n_vars + 1, 1] + 1:(length(n_ahead) * n_vars)

  # 'layout' reduces the size of the plotted text, which has to be taken into
  # account for the space that is required for the labels of the x-axis
  graphics::layout(layout_matrix)
  cex_lay <- graphics::par("cex")

  mar_ax_x <- graphics::par("mgp")[2] + label_width * cex_lay / csi + .3
  # The subplots must not become too small
  mar_ax_x <- min(mar_ax_x, in_height / csi / 3)

  # Space that is required for the axes, measured in fractions of the region
  # that is available for the subplots
  ax_width <- mar_ax_y * csi / in_width
  ax_height <- mar_ax_x * csi / in_height

  plot_width <- (1 - lab_size - ax_width) / length(n_ahead)
  plot_height <- (1 - lab_size - ax_height) / n_vars

  graphics::layout(layout_matrix,
                   widths = c(lab_size, ax_width + plot_width,
                              rep(plot_width, length(n_ahead) - 1)),
                   heights = c(lab_size, rep(plot_height, n_vars - 1),
                               ax_height + plot_height))

  # Forecast horizons above the columns of the plot
  for (j in 1:length(n_ahead)) {
    graphics::par(mar = c(0, ifelse(j == 1, mar_ax_y, mar_min), 0, mar_min))
    graphics::plot.new()
    graphics::text(0.5, 0.5, labels = paste0("h = ", n_ahead[j]), adj = 0.5)
  }

  # Variable names left of the rows of the plot
  for (i in 1:n_vars) {
    graphics::par(mar = c(ifelse(i == n_vars, mar_ax_x, mar_min), 0, mar_min, 0))
    graphics::plot.new()
    graphics::text(0.5, 0.5, labels = vars[i], adj = 0.5, srt = 90)
  }

  # The x-axis is not extended beyond 'xlim' so that the distance between the
  # outer error bars and the borders of a subplot is the same as the distance
  # between the error bars and the separators of the periods
  graphics::par(xaxs = "i")

  for (j in 1:length(n_ahead)) {
    for (i in 1:n_vars) {
      graphics::par(mar = c(ifelse(i == n_vars, mar_ax_x, mar_min),
                            ifelse(j == 1, mar_ax_y, mar_min),
                            mar_min, mar_min))

      graphics::plot(NA, xlim = xlim, ylim = ylim[[vars[i]]], axes = FALSE,
                     main = "", xlab = "", ylab = "", ...)

      # Separators between the periods
      if (length(periods) > 1) {
        graphics::abline(v = 1:(length(periods) - 1) + .5, lty = 3, col = "darkgrey")
      }

      if (zero_line) {
        graphics::abline(h = 0, col = "darkgrey")
      }

      for (m in 1:n_models) {
        temp <- errors[[m]]
        if (is.null(temp)) {
          # Forecast errors are not contained in the respective model
          next
        }
        temp <- temp[temp[, "variable"] == vars[i] & temp[, "h"] == n_ahead[j], ]
        if (nrow(temp) == 0) {
          # Variable or forecast horizon is not contained in the respective model
          next
        }

        pos <- match(temp[, "period"], periods) + offset[m]

        # The label of an error bar is the position of the model in 'x'
        draw_error_bars(pos, temp, cap = cap, labels = m,
                        col = col[m], lwd = lwd[m], pch = pch, cex = cex)
      }

      if (i == n_vars) {
        graphics::axis(1, at = 1:length(periods), labels = period_labels, las = 2)
      }
      if (j == 1) {
        graphics::axis(2)
      }
      graphics::box()
    }
  }

  if (!is.null(xlab)) {
    graphics::mtext(xlab, side = 1, outer = TRUE, line = .8,
                    cex = graphics::par("cex.lab") / graphics::par("cex"))
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = .5,
                    cex = graphics::par("cex.main") / graphics::par("cex"),
                    font = graphics::par("font.main"))
  }
}
