#' @include summary.bvec.R
#'
#' @export
#' @rdname summary.bvec
print.summary.bvec <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  k <- x$specifications$dims
  
  cat("\nModel:\n\n",
      paste("y ~ ", paste(dimnames(x$coefficients$means)[[2]],
                         collapse = " + "),
            sep = ""), "\n", sep = "")
  
  y_names <- dimnames(x$coefficients$means)[[1]]
  
  if (!is.null(x$coefficients)) {
    for (i in 1:k) {
      temp <- cbind(x$coefficients$means[i, ],
                    x$coefficients$sd[i, ],
                    x$coefficients$naivesd[i, ],
                    x$coefficients$tssd[i, ],
                    x$coefficients$q_lower[i, ],
                    x$coefficients$median[i, ],
                    x$coefficients$q_upper[i, ])
      
      dimnames(temp)[[2]] <- c("Mean", "SD", "Naive SD", "Time-series SD",
                               x$specifications$ci[1], "50%", x$specifications$ci[2])
      
      cat("\nVariable:", y_names[i], "\n\n")
      print.default(temp, quote = FALSE, right = TRUE, digits = digits)
    } 
  }
  
  if (!is.null(x$sigma)) {
    x_names <- NULL
    for (i in 1:k) {
      x_names <- c(x_names , paste(dimnames(x$sigma$means)[[1]][i], dimnames(x$sigma$means)[[1]], sep = "_"))
    }

    temp <- cbind(matrix(x$sigma$means),
                  matrix(x$sigma$sd),
                  matrix(x$sigma$naivesd),
                  matrix(x$sigma$tssd),
                  matrix(x$sigma$q_lower),
                  matrix(x$sigma$median),
                  matrix(x$sigma$q_upper))
    
    dimnames(temp) <- list(x_names, 
                           c("Mean", "SD", "Naive SD", "Time-series SD",
                             x$specifications$ci[1], "50%", x$specifications$ci[2]))
    
    if (k == 1) {
      cat("\nVariance:\n\n")
    } else {
      cat("\nVariance-covariance matrix:\n\n") 
    }
    print.default(temp, quote = FALSE, right = TRUE, digits = digits)
  }
  
  cat("\n")
  invisible(x)
}
