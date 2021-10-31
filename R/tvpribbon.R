
# Use for the extraction of time-varying parameters in plot.bvar()
.tvpribbon <- function(x, var, ymin, ymax) {
  draws <- ts(cbind(t(matrix(unlist(lapply(x, function(x, var, ymin, ymax) {quantile(x[, var], probs = c(ymin, .5, ymax))},
                                           var = var, ymin = ymin, ymax = ymax)), 3)), 0))
}