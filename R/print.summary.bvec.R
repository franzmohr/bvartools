#' @include summary.bvec.R
#'
#' @export
#' @rdname summary.bvec
print.summary.bvec <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  cat("\nModel:\n\n",
      paste("y ~ ", paste(dimnames(x$coefficients$means)[[2]],
                         collapse = " + "),
            sep = ""), "\n", sep = "")
  
  k <- x$specifications$dims["K"]
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
      
      dim_names_1 <- dimnames(x$coefficients$means)[[2]]
      dim_names_2 <- c("Mean", "SD", "Naive SD", "Time-series SD",
                       x$specifications$ci[1], "50%", x$specifications$ci[2])
      
      if ("lambda" %in% names(x$coefficients)) {
        temp <- cbind(temp, x$coefficients$lambda[i, ])
        dim_names_2 <- c(dim_names_2, "Incl. prob.")
      }
      
      dimnames(temp)[[1]] <- dim_names_1
      dimnames(temp)[[2]] <- dim_names_2
      
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
    
    dim_names_1 <- x_names
    dim_names_2 <- c("Mean", "SD", "Naive SD", "Time-series SD",
                     x$specifications$ci[1], "50%", x$specifications$ci[2])
    
    
    if ("lambda" %in% names(x$sigma)) {
      temp <- cbind(temp, matrix(x$sigma$lambda))
      dim_names_2 <- c(dim_names_2, "Incl. prob.")
    }
    
    dimnames(temp) <- list(dim_names_1, 
                           dim_names_2)
    
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
