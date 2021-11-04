
# Use for the extraction of time-varying parameters in plot.bvar()
.tvpribbon <- function(x, var, ymin, ymax) {
  draws <- stats::ts(cbind(t(matrix(unlist(lapply(x, function(x, var, ymin, ymax) {stats::quantile(x[, var], probs = c(ymin, .5, ymax))},
                                           var = var, ymin = ymin, ymax = ymax)), 3)), 0))
}