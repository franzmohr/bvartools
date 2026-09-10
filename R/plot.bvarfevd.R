#' Plotting Forecast Error Variance Decompositions of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarfevd".
#' 
#' @param x an object of class "bvarfevd", usually, a result of a call to \code{\link{fevd}}.
#' @param max_groups integer. Maximum number of variables shown in the plot. The
#' \code{max_groups - 1} variables with the largest contributions across the whole horizon
#' are kept and the contributions of the remaining variables are added up in a further bar
#' segment named \code{"Other"}, so that the legend does not become too large. Default is
#' \code{NULL}, so that a segment is shown for every variable of the decomposition.
#' @param ... further graphical parameters.

#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- add_posterior_coefficients(model)
#' 
#' # Obtain FEVD
#' vd <- fevd(object, response = "cons")
#' 
#' # Plot
#' plot(vd)
#' 
#' @export
#' @rdname fevd
plot.bvarfevd <- function(x, max_groups = NULL, ...) {
  # Only save and restore 'mar'. Restoring all parameters would also reset
  # 'mfg' and thus overwrite the current panel of a user-defined layout.
  orig_par <- graphics::par(mar = graphics::par("mar"))
  on.exit(graphics::par(orig_par))

  periods <- stats::time(x)

  # Pool the smallest contributions into one segment, if asked for. Done on a
  # plain matrix, since the time index is already taken and subsetting a "ts"
  # drops the column names the legend is built from.
  max_groups <- .check_max_groups(max_groups)
  shares <- unclass(x)
  attr(shares, "tsp") <- NULL
  if (!is.null(max_groups) && max_groups < ncol(shares)) {
    shares <- .limit_fevd_groups(shares, max_groups)
  }

  legend_names <- colnames(shares)

  char_width <- graphics::par("cin")[1] # Width of a character in inches
  line_height <- graphics::par("csi")   # Height of a line of text in inches

  # Space needed by the legend (labels, fill boxes and padding) in inches and
  # the space that may be spent on it without squashing the bars
  legend_width <- (max(nchar(legend_names)) + 2.5) * char_width + .3
  max_width <- .4 * graphics::par("fin")[1]

  # Shrink the legend text if the labels do not fit into the available space
  legend_cex <- 1
  if (legend_width > max_width) {
    legend_cex <- max(max_width / legend_width, .5)
    legend_width <- min(legend_width * legend_cex, max_width)
  }

  graphics::par(mar = c(5.1, 4.1, 4.1, legend_width / line_height))

  # Defaults that may be overridden via '...'
  dots <- list(...)
  args <- list(ylab = "Percentage", xlab = "Period", names.arg = periods)
  args <- args[!(names(args) %in% names(dots))]

  do.call(graphics::barplot, c(list(t(shares)), args, dots))

  # Place the legend in the right margin, vertically centred on the plot region
  usr <- graphics::par("usr")
  graphics::legend(x = usr[2] + .02 * (usr[2] - usr[1]), y = mean(usr[3:4]),
                   xjust = 0, yjust = .5, xpd = TRUE, cex = legend_cex,
                   legend = legend_names, fill = grDevices::gray.colors(ncol(shares)))
}
