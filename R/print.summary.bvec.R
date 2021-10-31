#' @include summary.bvec.R
#'
#' @export
#' @rdname summary.bvec
print.summary.bvec <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  # Title
  title_text <- "\nBayesian "
  tvp <- any(unlist(x[["specifications"]][["tvp"]])[c("A0", "Pi", "Pi_x", "Pi_d", "Gamma", "Upsilon", "C")])
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (x[["specifications"]][["tvp"]][["Sigma"]]) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["specifications"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VEC model")
  p_text <- paste0("p = ", x[["specifications"]][["lags"]][["p"]])
  s_text <- NULL
  if (x[["specifications"]][["dims"]][["M"]] > 0) {
    s_text <- paste0("s = ", x[["specifications"]][["lags"]][["s"]])
  }
  if (any(!is.null(c(p_text, s_text)))) {
    lag_text <- paste0(c(p_text, s_text), collapse = " and ")
  } else {
    lag_text <- NULL
  }
  title_text <- paste0(c(title_text, lag_text), collapse = " with ")
  
  cat(title_text, "\n")
  
  cat("\nModel:\n\n",
      paste("d.y ~ ", paste(dimnames(x$coefficients$means)[[2]],
                         collapse = " + "),
            sep = ""), "\n", sep = "")
  
  
  if (!is.null(x[["specifications"]][["period"]])) {
    cat("\nPeriod:", x[["specifications"]][["period"]], "\n")
  }
  
  k <- x[["specifications"]][["dims"]][["K"]]
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
