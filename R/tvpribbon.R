
# Use for the extraction of time-varying parameters in plot.bvar()
.tvpribbon <- function(x, var, ymin, ymax, show_zero_y) {
  draws <- stats::ts(t(matrix(unlist(lapply(x, function(x, var, ymin, ymax) {stats::quantile(x[, var], probs = c(ymin, .5, ymax))},
                                            var = var, ymin = ymin, ymax = ymax)), 3)))
  
  if (show_zero_y) {
    draws <- cbind(draws, 0)
  }
  return(draws)
}