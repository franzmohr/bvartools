#' Plotting Forecast Error Variance Decompositions of Bayesian Vector Autoregression
#' 
#' A plot function for objects of class "bvarfevd".
#' 
#' @param x an object of class "bvarfevd", usually, a result of a call to \code{\link{fevd}}.
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
plot.bvarfevd <- function(x, ...) {
  # Only save and restore 'mar'. Restoring all parameters would also reset
  # 'mfg' and thus overwrite the current panel of a user-defined layout.
  orig_par <- graphics::par(mar = graphics::par("mar"))
  on.exit(graphics::par(orig_par))

  legend_names <- dimnames(x)[[2]]

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
  args <- list(ylab = "Percentage", xlab = "Period", names.arg = stats::time(x))
  args <- args[!(names(args) %in% names(dots))]

  do.call(graphics::barplot, c(list(t(x)), args, dots))

  # Place the legend in the right margin, vertically centred on the plot region
  usr <- graphics::par("usr")
  graphics::legend(x = usr[2] + .02 * (usr[2] - usr[1]), y = mean(usr[3:4]),
                   xjust = 0, yjust = .5, xpd = TRUE, cex = legend_cex,
                   legend = legend_names, fill = grDevices::gray.colors(NCOL(x)))
}