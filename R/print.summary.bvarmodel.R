#' @include summary.bvarmodel.R
#'
#' @export
#' @rdname summary.bvarmodel
print.summary.bvarmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  k <- x[["model"]][["k"]]
  
  # Title
  title_text <- "\nBayesian " 
  tvp <- x[["model"]][["tvp"]]
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (x[["model"]][["error"]] %in% c("sv", "sv+covar")) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  if (k == 1) {
    title_text <- paste0(title_text, "AR model")
  } else {
    title_text <- paste0(title_text, "VAR model") 
  }
  p_text <- paste0("p = ", x[["model"]][["p"]])
  s_text <- NULL
  if (x[["model"]][["m"]] > 0) {
    s_text <- paste0("s = ", x[["model"]][["s"]])
  }
  if (any(!is.null(c(p_text, s_text)))) {
    lag_text <- paste0(c(p_text, s_text), collapse = " and ")
  } else {
    lag_text <- NULL
  }
  title_text <- paste0(c(title_text, lag_text), collapse = " with ")
  
  cat(title_text, "\n")
  
  use_incl <- x[["model"]][["varsel"]] %in% c("ssvs", "bvs")
  if (use_incl) {
    if (x[["model"]][["varsel"]] == "bvs") {
      varsel_algo <- "Bayesian variable selection (Korobilis, 2013)"
    }
    if (x[["model"]][["varsel"]] == "ssvs") {
      varsel_algo <- "Stochastic search variable selection (George et al., 2008)"
    }
    cat("\nVariable selection algorithm: ", varsel_algo, "\n", sep = "")
  }
  
  # Model
  
  if (is.null(x[["a"]][["means"]])) {
    regressors <- "nothing"
    use_a <- FALSE
  } else {
    regressors <- paste(dimnames(x[["a"]][["means"]])[[1]], collapse = ", ")
    use_a <- TRUE
  }
  
  if (k == 1) {
    cat(paste("\nEndogenous variable: ", regressors, sep = ""), "\n", sep = "")
  } else {
    cat(paste("\nEndogenous variables: ", regressors, sep = ""), "\n", sep = "") 
  }
  
  if (!is.null(x[["model"]][["period"]])) {
    cat("\nPeriod:", x[["model"]][["period"]], "\n")
  }
  
  
  y_names <- x[["model"]][["endogen"]]
  
  # Coefficients per endogenous variable
  
  if (use_a) {
    for (i in 1:k) {
      temp <- cbind(x[["a"]][["means"]][i, ],
                    x[["a"]][["sd"]][i, ],
                    x[["a"]][["naivesd"]][i, ],
                    x[["a"]][["tssd"]][i, ],
                    x[["a"]][["q_lower"]][i, ],
                    x[["a"]][["median"]][i, ],
                    x[["a"]][["q_upper"]][i, ])
      
      dim_names_1 <- dimnames(x[["a"]][["means"]])[[2]]
      dim_names_2 <- c("Mean", "SD", "Naive SD", "Time-series SD",
                       x[["model"]][["ci"]][1], "50%", x[["model"]][["ci"]][2])
      
      if ("lambda" %in% names(x[["a"]])) {
        temp <- cbind(temp, x[["a"]][["lambda"]][i, ])
        dim_names_2 <- c(dim_names_2, "Incl. prob.")
      }
      
      dimnames(temp)[[1]] <- dim_names_1
      dimnames(temp)[[2]] <- dim_names_2
      
      if (k > 1) {
        cat("\nVariable:", y_names[i], "\n\n") 
      }
      temp <- as.data.frame(temp)
      temp[, "Signif."] <- sign(temp[, x[["model"]][["ci"]][1]]) == sign(temp[, x[["model"]][["ci"]][2]])
      temp[, "Signif."] <- ifelse(temp[, "Signif."], "*", "")
      names(temp)[length(names(temp))] <- ""
      print(temp, ...)
    } 
  } else {
    cat("\n\nNo regressors.\n\n")
  }
  
  # Error covariance matrix
  
  if (!is.null(x[["sigma"]])) {
    x_names <- NULL
    for (i in 1:k) {
      x_names <- c(x_names , paste(dimnames(x[["sigma"]][["means"]])[[1]][i], dimnames(x[["sigma"]][["means"]])[[1]], sep = "_"))
    }
    
    temp <- cbind(matrix(x[["sigma"]][["means"]]),
                  matrix(x[["sigma"]][["sd"]]),
                  matrix(x[["sigma"]][["naivesd"]]),
                  matrix(x[["sigma"]][["tssd"]]),
                  matrix(x[["sigma"]][["q_lower"]]),
                  matrix(x[["sigma"]][["median"]]),
                  matrix(x[["sigma"]][["q_upper"]]))
    
    dim_names_1 <- x_names
    dim_names_2 <- c("Mean", "SD", "Naive SD", "Time-series SD",
                     x[["model"]][["ci"]][1], "50%", x[["model"]][["ci"]][2])
    
    if ("lambda" %in% names(x[["sigma"]])) {
      temp <- cbind(temp, matrix(x[["sigma"]][["lambda"]]))
      dim_names_2 <- c(dim_names_2, "Incl. prob.")
    }
    
    dimnames(temp) <- list(dim_names_1, 
                           dim_names_2)
    
    temp <- temp[which(lower.tri(matrix(1:(k * k), k), diag = TRUE)), ]
    
    if (k == 1) {
      cat("\nVariance:\n\n")
    } else {
      cat("\nVariance-covariance matrix:\n\n") 
    }
    temp <- as.data.frame(temp)
    temp[, "Signif."] <- sign(temp[, x[["model"]][["ci"]][1]]) == sign(temp[, x[["model"]][["ci"]][2]])
    temp[, "Signif."] <- ifelse(temp[, "Signif."], "*", "")
    names(temp)[length(names(temp))] <- ""
    print(temp, ...)
  }
  
  cat("\n")
  invisible(x)
}
