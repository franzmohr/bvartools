#' Plotting Factors from Dynamic Factor Models
#' 
#' A plot function for objects of class \code{"dfm"}.
#' 
#' @param x an object of class \code{"dfm"}, usually, a result of a call to \code{\link{dfm}}.
#' @param ci interval used to calculate credible bands.
#' @param ... further graphical parameters.
#' 
#' @examples
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(v_i = .01),
#'                     sigma_u = list(shape = 5, rate = 4),
#'                     a = list(v_i = .01),
#'                     sigma_v = list(shape = 5, rate = 4))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Plot factors
#' plot(object)
#' 
#' @export
#' @rdname dfm
plot.dfm <- function(x, ci = 0.95, ...) {
  
  m <- x[["specifications"]][["dims"]]["M"]
  n <- x[["specifications"]][["dims"]]["N"]
  tt <- ncol(x[["factor"]]) / n
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(x[["factor"]], 2, stats::quantile, probs = c(ci_low, .5, ci_high))
  median <- NULL
  q_low <- NULL
  q_high <- NULL
  for (i in 1:n) {
    q_low <- cbind(q_low, matrix(temp[1, i + n * 0:(tt - 1)], tt))
    median <- cbind(median, matrix(temp[2, i + n * 0:(tt - 1)], tt))
    q_high <- cbind(q_high, matrix(temp[3, i + n * 0:(tt - 1)], tt)) 
  }
  
  var_names <- paste("Factor", 1:n)
  dimnames(q_low) <- list(NULL, var_names)
  dimnames(median) <- list(NULL, var_names)
  dimnames(q_high) <- list(NULL, var_names)
  
  graphics::par(mfcol = c(length(var_names), 1))
  for (i in var_names) {
    temp <- cbind(q_low[, i], median[, i], q_high[, i])
    temp <- stats::ts(temp)
    stats::tsp(temp) <- stats::tsp(x$x)
    stats::plot.ts(temp, plot.type = "single", lty = c(2, 1, 2), main = i, ylab = "")
  }
  graphics::par(mfcol = c(1, 1))
}
