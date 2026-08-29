#' Plotting Forecast Errors
#' 
#' A plot function for objects of class 'fcsterrorlist'.
#' 
#' @param x an object of class 'fcsterrorlist', usually, a result of a call to \code{\link{forecast_errors}}.
#' @param ... further graphical parameters.
#' 
#' @export
plot.fcsterrorlist <- function(x, ...) {
  
  # 'layout' is called below, so all parameters have to be restored on exit
  orig_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(orig_par))

  n_models <- length(x)
  vars <- unique(unlist(lapply(x, function(x) {names(x)})))
  n_vars <- length(vars)
  n_ahead <- unique(unlist(lapply(x, function(x) {dim(x[[1]])[2]})))
  
  lab_size <- .05
  
  layout_matrix <- matrix(0, nrow = n_ahead + 1, ncol = length(vars) + 1)
  layout_matrix[1, -1] <- 1:length(vars)
  layout_matrix[-1, 1] <- length(vars) + 1:n_ahead
  layout_matrix[-1, -1] <- layout_matrix[n_ahead + 1, 1] + 1:(length(vars) * n_ahead)
  graphics::layout(layout_matrix,
                   widths = c(lab_size, rep((1 - lab_size) / n_vars, n_vars)),
                   heights = c(lab_size, rep((1 - lab_size) / n_ahead, n_ahead)))
  
  graphics::par(mar = c(0, 0, 0, 0))
  for (j in vars) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  
  graphics::par(mar = c(3, 0, 0, 0))
  for (j in 1:n_ahead) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = paste0("h = ", j), adj = 0.5)
  }
  
  graphics::par(mar = c(3, 1.5, 0, 1))
  
  for (var in vars) {
    for (h in 1:n_ahead) {
      temp <- list()
      for (mod in 1:n_models) {
        temp[[mod]] <- abs(stats::na.omit(c(x[[mod]][[var]][, h, ])))
      }
      graphics::boxplot(temp, ...)
    }
  }
  
}


