#' Plotting Selection Criteria
#' 
#' A plot function for objects of class 'selcritlist'.
#' 
#' @param x an object of class 'selcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are the in-sample criteria \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default),
#' \code{"HQ"} and the out-of-sample statistics \code{"FE"}, \code{"AFE"} and
#' \code{"RSFE"}.
#' @param ... further graphical parameters.
#'
#' @details For in-sample criteria each model is represented by a horizontal error
#' bar, which covers the credible band of the criterion and which is marked by the
#' median and the mean of its posterior draws. The criterion is measured on the
#' x-axis and the models are arranged along the y-axis, beginning with the first
#' model of \code{x} at the top. This keeps the plot readable, if \code{x} contains
#' many models. The position of a model in \code{x} is added at the upper end of its
#' error bar.
#'
#' The specifications of the models are used as labels of the y-axis, where only the
#' specifications that differ across the models are used, as in
#' \code{\link{print.modellist}}. The lag order is always used. The left margin is
#' adjusted to the width of the labels.
#'
#' No title is added by default. Space above the plot is only used, if argument
#' \code{main} is specified.
#'
#' For out-of-sample statistics each row of the plot corresponds to an endogenous
#' variable and the x-axis to the forecast horizon. Within a forecast horizon each
#' model is represented by an error bar, which covers the credible band of the
#' respective statistic and which is marked by the median and the mean of its
#' posterior draws. The position of a model in \code{x} is added above its error
#' bar. Since forecast errors can be positive and negative, a horizontal reference
#' line is added at zero for \code{criterion = "FE"}.
#'
#' For both types of criteria arguments \code{col}, \code{pch}, \code{cex},
#' \code{lwd}, \code{main} and \code{xlab} can be used to change the appearance of
#' the plot, where \code{col} and \code{lwd} are recycled over the models in
#' \code{x}.
#'
#' @export
plot.selcritlist <- function(x, criterion = "BIC", ...) {

  in_sample <- c("LL", "AIC", "BIC", "HQ")
  out_of_sample <- c("FE", "AFE", "RSFE")

  avail_statistics <- unique(unlist(lapply(x, function(y) {names(y)})))
  avail_statistics <- avail_statistics[avail_statistics %in% c(in_sample, out_of_sample)]

  if (length(avail_statistics) == 0) {
    stop("Function only supports ", paste0(c(in_sample, out_of_sample), collapse = ", "),
         ". None is contained in argument 'x'.")
  }

  if (!criterion %in% avail_statistics) {
    stop("Criterion '", criterion, "' is not contained in argument 'x'. Available criteria are ",
         paste0("'", avail_statistics, "'", collapse = ", "), ".")
  }

  # Out-of-sample statistics are plotted per variable and forecast horizon
  if (criterion %in% out_of_sample) {
    plot_selcrit_forecast_errors(x, criterion = criterion, ...)
    return(invisible(x))
  }

  n_models <- length(x)

  res <- lapply(x, function(y, criterion) {y[[criterion]]}, criterion = criterion)
  res <- do.call("rbind", res)[, c("mean", "median", "qlower", "qupper")]

  # Model specifications for x-axis
  specs <- lapply(x, get_model_specifications)

  if (length(unique(lapply(specs, names))) > 1) {
    stop("The models in argument 'x' do not contain the same specifications, so they cannot be compared.")
  }
  specs <- do.call("rbind", specs)

  # Models of class 'bvecmodel' distinguish between restricted and unrestricted
  # deterministic terms, which are combined for the labels
  if (all(c("n_unrestricted", "n_restricted") %in% names(specs))) {
    specs[["n_unrestricted"]] <- specs[["n_unrestricted"]] + specs[["n_restricted"]]
    specs[["n_restricted"]] <- NULL
    names(specs)[names(specs) == "n_unrestricted"] <- "n"
  }

  # Only variables, which differ across the models, are used, because identical
  # values do not help to distinguish the models. The lag order is always used
  keep <- unlist(lapply(specs, function(y) {any(y != y[1], na.rm = TRUE)}))
  if ("p" %in% names(keep)) {
    keep["p"] <- TRUE
  }
  specs <- specs[, keep, drop = FALSE]

  # Each element of 'model_specs' contains the lines of the label of a model
  model_specs <- list()
  for (i in 1:n_models) {
    temp <- NULL
    for (j in names(specs)) {
      if (j == "type") {
        temp <- c(temp, as.character(specs[i, j]))
      } else if (j == "varsel") {
        temp <- c(temp, paste0("(", specs[i, j], ")"))
      } else {
        temp <- c(temp, paste0(ifelse(j == "rank", "r", j), "=", specs[i, j]))
      }
    }
    model_specs[[i]] <- temp
  }

  # The model descriptions are used as labels of the y-axis, so all
  # specifications of a model are combined in a single line
  labels <- unlist(lapply(model_specs, paste0, collapse = ", "))

  dots <- list(...)

  # The axes and the title are drawn by the function, so the respective
  # arguments are dropped
  xlab <- dots[["xlab"]]
  dots[["xlab"]] <- NULL
  ylab <- dots[["ylab"]]
  dots[["ylab"]] <- NULL
  dots[["yaxt"]] <- NULL
  main <- dots[["main"]]
  dots[["main"]] <- NULL
  # The criterion is measured on the x-axis
  if (is.null(xlab)) {xlab <- criterion}
  if (is.null(ylab)) {ylab <- ""}
  # Arguments of the error bars
  col <- dots[["col"]]
  dots[["col"]] <- NULL
  lwd <- dots[["lwd"]]
  dots[["lwd"]] <- NULL
  pch <- dots[["pch"]]
  dots[["pch"]] <- NULL
  cex_pt <- dots[["cex"]]
  dots[["cex"]] <- NULL
  if (is.null(col)) {col <- "black"}
  if (is.null(lwd)) {lwd <- 1}
  if (is.null(pch)) {pch <- 20}
  if (is.null(cex_pt)) {cex_pt <- 1}
  col <- rep_len(col, n_models)
  lwd <- rep_len(lwd, n_models)

  cex_axis <- ifelse(is.null(dots[["cex.axis"]]), graphics::par("cex.axis"), dots[["cex.axis"]])
  cex_lab <- ifelse(is.null(dots[["cex.lab"]]), graphics::par("cex.lab"), dots[["cex.lab"]])
  cex_main <- ifelse(is.null(dots[["cex.main"]]), graphics::par("cex.main"), dots[["cex.main"]])

  orig_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = orig_mar))
  mar <- orig_mar

  # Space that is required for the labels of the axes, measured in lines of text
  csi <- graphics::par("csi")
  line_ticks <- graphics::par("mgp")[2]
  lab_extent_y <- line_ticks +
    max(graphics::strwidth(labels, units = "inches", cex = cex_axis)) / csi
  lab_extent_x <- line_ticks + cex_axis

  # The margins only contain the labels of the axes. Space above the plot is
  # only used, if a title is added
  line_ylab <- lab_extent_y + .3
  line_xlab <- lab_extent_x + .5
  mar[1] <- line_xlab + ifelse(nzchar(xlab), cex_lab + .3, 0)
  mar[2] <- line_ylab + ifelse(nzchar(ylab), cex_lab + .3, 0)
  # The plot region must not become too small
  mar[2] <- min(mar[2], graphics::par("din")[1] / csi / 2)
  mar[3] <- ifelse(is.null(main), .5, 1 + cex_main)
  mar[4] <- .5
  graphics::par(mar = mar)

  # Common limits of the axes. The first model is plotted at the top
  ylim <- c(n_models, 0) + c(.5, .5)
  xlim <- range(unlist(res)[is.finite(unlist(res))])
  if (diff(xlim) == 0) {
    xlim <- xlim + c(-.5, .5)
  }
  # Space for the model numbers at the upper end of the error bars
  xlim <- xlim + c(0, .1 * diff(xlim))

  # Half of the width of the lines, which mark the bounds and the median
  cap <- min(.2, n_models / 10)

  # Defaults that may be overridden via '...'. The y-axis is not extended beyond
  # 'ylim' so that the distance between the outer error bars and the borders of
  # the plot is the same as the distance between the error bars
  args <- list(yaxs = "i")
  args <- args[!(names(args) %in% names(dots))]

  do.call(graphics::plot, c(list(NA), list(xlim = xlim, ylim = ylim, axes = FALSE,
                                          main = "", xlab = "", ylab = ""), args, dots))

  draw_error_bars(1:n_models, res, cap = cap, labels = 1:n_models, horizontal = TRUE,
                  col = col, lwd = lwd, pch = pch, cex = cex_pt)

  graphics::axis(1)
  graphics::axis(2, at = 1:n_models, labels = labels, las = 1, cex.axis = cex_axis)
  graphics::box()

  if (nzchar(xlab)) {
    graphics::mtext(xlab, side = 1, line = line_xlab, cex = cex_lab)
  }
  if (nzchar(ylab)) {
    graphics::mtext(ylab, side = 2, line = line_ylab, cex = cex_lab)
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, line = .5, cex = cex_main,
                    font = graphics::par("font.main"))
  }

  return(invisible(x))
}
