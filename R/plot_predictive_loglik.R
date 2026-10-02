#' Plotting the Log Predictive Likelihood per Period
#'
#' Plots the log predictive density of each evaluated period for one or more
#' Bayesian models, either period by period or as a running sum.
#'
#' @param x an object of class 'selcrit' or 'selcritlist' whose entries contain
#' the log predictive likelihood, usually a result of a call to
#' \code{\link{selection_criteria}} on expanding windows that carry the draws of
#' \code{\link{add_predictive_loglik}}.
#' @param baseline the model every other model is compared with, given as its
#' position in \code{x} or as its name. If \code{NULL} (default), nothing is
#' subtracted and the log predictive densities themselves are plotted.
#' @param cumulative logical. Should the running sum over the periods be plotted
#' (\code{TRUE}, the default) rather than the value of each period?
#' @param col a vector of colours, which is recycled over the models in \code{x}.
#' @param lwd a vector of line widths, which is recycled over the models in
#' \code{x}.
#' @param lty a vector of line types, which is recycled over the models in
#' \code{x}.
#' @param legend logical. Should a legend of the models be added?
#' @param legend_position the position of the legend, passed on to
#' \code{\link[graphics]{legend}}. The default is \code{"topleft"}, which a
#' cumulative plot of differences often leaves free; where the lines rise into
#' it, \code{"bottomleft"} or \code{"right"} usually does not.
#' @param main the title of the plot. If \code{NULL} (default), no title is added.
#' @param xlab the label of the x-axis.
#' @param ylab the label of the y-axis. If \code{NULL} (default), a label is
#' chosen from \code{baseline} and \code{cumulative}.
#' @param ... further graphical parameters, which are passed on to
#' \code{\link[graphics]{plot}}.
#'
#' @details The log predictive likelihood of a model is the sum of the log
#' predictive densities of the periods it was evaluated for, and
#' \code{\link{selection_criteria}} reports that sum with a band. This function
#' plots the terms of the sum, which answers a question the sum cannot: whether a
#' difference between two models is a steady advantage or the consequence of a
#' few periods. A single quarter can decide a comparison of this kind, and with
#' \code{cumulative = TRUE} such a quarter is a step in an otherwise flat line
#' rather than a number hidden in a total.
#'
#' \strong{With a \code{baseline} the plot shows differences}, the cumulative sum
#' of \eqn{lpd_{t} - lpd_{t}^{baseline}} over the evaluated periods, and a
#' horizontal reference line is drawn at zero. A line that ends above zero
#' belongs to a model that predicted better than the baseline over the whole
#' evaluation sample, and the slope at any point is how much better it was
#' predicting then. Differences are what the criterion is read in, so a baseline
#' is usually what is wanted as soon as there is more than one model; without one
#' the lines are dominated by the periods that were hard for every model alike.
#'
#' The periods come from attribute \code{"terms"} of the \code{LPL} entry of each
#' model, which also holds the numerical standard error of each period. Models are
#' plotted over the periods they share with the baseline, and a model without an
#' \code{LPL} entry is left out with a warning.
#'
#' @return \code{x}, invisibly. The function is called for its side effect, the
#' plot.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' # Two specifications to compare
#' models <- create_bvarmodel(data = us_macrodata, p = 1:2,
#'                            deterministic = "const",
#'                            iterations = 100, burnin = 50)
#' models <- add_priors(models,
#'                      coef = list(v_i = 1, v_i_det = 1 / 10),
#'                      sigma = list(df = 3, scale = 1))
#'
#' # One window per evaluated period
#' windows <- lapply(models, use_expanding_window, start = c(1972, 1))
#'
#' criteria <- lapply(windows, function(w) {
#'   w <- add_initial_values(w)
#'   w <- add_posterior_coefficients(w)
#'   w <- add_predictive_loglik(w)
#'   selection_criteria(w)
#' })
#' class(criteria) <- append("selcritlist", class(criteria))
#'
#' plot_predictive_loglik(criteria, baseline = 1)
#'
#' @seealso \code{\link{selection_criteria}} for the sum and its band,
#' \code{\link{add_predictive_loglik}} for the draws it is computed from, and
#' \code{\link{plot_forecast_errors_by_period}} for the same question asked of
#' forecast errors.
#'
#' @export
plot_predictive_loglik <- function(x, baseline = NULL, cumulative = TRUE,
                                   col = "black", lwd = 1, lty = 1,
                                   legend = TRUE, legend_position = "topleft",
                                   main = NULL, xlab = "Period",
                                   ylab = NULL, ...) {

  if (!is.logical(cumulative) || length(cumulative) != 1 || is.na(cumulative)) {
    stop("Argument 'cumulative' must be TRUE or FALSE.")
  }

  object <- x

  if (!any(c("selcrit", "selcritlist") %in% class(x))) {
    stop("Argument 'x' must be an object of class 'selcrit' or 'selcritlist'. ",
         "Use selection_criteria() on expanding windows that carry the draws of ",
         "add_predictive_loglik().")
  }

  # A single object of class 'selcrit' is one model
  if (!"selcritlist" %in% class(x)) {
    x <- list(x)
  }

  # The periods and their log predictive densities, one data frame per model.
  # Models without an LPL entry are dropped rather than silently plotted as a
  # gap, since 'baseline' counts positions in 'x'.
  terms_of <- function(entry) {
    if (is.null(entry) || is.null(entry[["LPL"]])) {
      return(NULL)
    }
    attr(entry[["LPL"]], "terms")
  }
  terms <- lapply(x, terms_of)

  labels <- names(x)
  if (is.null(labels)) {
    labels <- paste("Model", seq_along(x))
  }

  missing <- vapply(terms, is.null, logical(1))
  if (all(missing)) {
    stop("No model in argument 'x' contains the log predictive likelihood. ",
         "add_predictive_loglik() has to be called on the windows before ",
         "selection_criteria().")
  }
  if (any(missing)) {
    warning("The log predictive likelihood is missing for ",
            paste0("'", labels[missing], "'", collapse = ", "),
            " and those models are not plotted.")
  }

  # The baseline, as a position in 'x'
  base_position <- NULL
  if (!is.null(baseline)) {
    if (length(baseline) != 1) {
      stop("Argument 'baseline' must be a single position or name.")
    }
    base_position <- if (is.character(baseline)) {
      match(baseline, labels)
    } else {
      as.integer(baseline)
    }
    if (is.na(base_position) || base_position < 1 || base_position > length(x)) {
      stop("Argument 'baseline' does not refer to a model in 'x'.")
    }
    if (missing[base_position]) {
      stop("The model given as 'baseline' does not contain the log predictive ",
           "likelihood.")
    }
  }

  n_models <- length(x)
  col <- rep_len(col, n_models)
  lwd <- rep_len(lwd, n_models)
  lty <- rep_len(lty, n_models)

  # One series per model: the period and the value plotted, which is the log
  # predictive density, the difference to the baseline, or either of them
  # accumulated over the periods.
  series <- vector("list", n_models)
  for (i in seq_len(n_models)) {
    if (missing[i]) {
      next
    }
    d <- terms[[i]][, c("period", "lpd")]
    d <- d[order(d[["period"]]), ]

    if (!is.null(base_position)) {
      b <- terms[[base_position]][, c("period", "lpd")]
      d <- merge(d, b, by = "period", suffixes = c("", ".base"))
      d <- d[order(d[["period"]]), ]
      d[["value"]] <- d[["lpd"]] - d[["lpd.base"]]
    } else {
      d[["value"]] <- d[["lpd"]]
    }

    if (nrow(d) == 0) {
      next
    }
    if (cumulative) {
      d[["value"]] <- cumsum(d[["value"]])
    }
    series[[i]] <- d[, c("period", "value")]
  }

  drawn <- which(!vapply(series, is.null, logical(1)))
  if (length(drawn) == 0) {
    stop("No model shares an evaluated period with the baseline.")
  }

  xlim <- range(unlist(lapply(series[drawn], function(d) d[["period"]])))
  ylim <- range(unlist(lapply(series[drawn], function(d) d[["value"]])))
  # The zero line is part of the message of a difference, so it is always inside
  # the window when there is a baseline.
  if (!is.null(base_position)) {
    ylim <- range(c(ylim, 0))
  }

  if (is.null(ylab)) {
    ylab <- if (is.null(base_position)) {
      if (cumulative) "Cumulative log predictive density" else "Log predictive density"
    } else {
      if (cumulative) {
        paste0("Cumulative difference to ", labels[base_position])
      } else {
        paste0("Difference to ", labels[base_position])
      }
    }
  }

  graphics::plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                 main = main, ...)
  if (!is.null(base_position)) {
    graphics::abline(h = 0, col = "grey60")
  }
  for (i in drawn) {
    graphics::lines(series[[i]][["period"]], series[[i]][["value"]],
                    col = col[i], lwd = lwd[i], lty = lty[i])
  }
  if (isTRUE(legend) && length(drawn) > 1) {
    graphics::legend(legend_position, legend = labels[drawn], col = col[drawn],
                     lwd = lwd[drawn], lty = lty[drawn], bty = "n")
  }

  return(invisible(object))
}
