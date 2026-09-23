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
  if (x[["model"]][["error"]] == "ald") {
    title_text <- paste0(title_text, "Quantile-")
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
  # The quantile is part of what the model is, not of how it was fitted, so it
  # is reported beside the lag orders rather than left to the specification.
  q_text <- NULL
  if (!is.null(x[["model"]][["quantile"]])) {
    q_text <- paste0("q = ", x[["model"]][["quantile"]])
  }
  if (any(!is.null(c(p_text, s_text, q_text)))) {
    lag_text <- paste0(c(p_text, s_text, q_text), collapse = " and ")
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

  # A sign restricted identification is set valued and its sample is not the
  # posterior sample, so both the restrictions and how many draws survived them
  # are reported. An acceptance rate is not a diagnostic of the sampler: a low
  # one says that the model rarely produces the pattern that was asked of it,
  # which is a finding about the restrictions.
  sign_restrictions <- x[["model"]][["sign_restrictions"]]
  if (!is.null(sign_restrictions)) {
    varnames <- x[["model"]][["endogen"]]
    printed <- sign_restrictions[["restrictions"]]
    printed[["impulse"]] <- varnames[printed[["impulse"]]]
    printed[["response"]] <- varnames[printed[["response"]]]
    printed[["sign"]] <- c("-", "0", "+")[printed[["sign"]] + 2]

    arw <- x[["model"]][["sign_zero_restrictions"]]
    cat("\nSign", if (is.null(arw)) NULL else "and zero", "restrictions:\n\n")
    print(printed, row.names = FALSE)

    if (is.null(arw)) {
      accepted <- sign_restrictions[["accepted"]]
      draws <- sign_restrictions[["draws"]]
      if (!is.null(accepted) && !is.null(draws) && draws > 0) {
        cat("\nIdentified draws: ", accepted, " of ", draws, " (",
            format(round(100 * accepted / draws, 1), nsmall = 1), "%)\n", sep = "")
      }
    } else {
      # The draws of an importance sample are resampled, so counting the
      # identified ones would report the resample rather than what went into
      # it. What says how much independent information is there is the
      # effective sample size, which the algorithm's own authors ask for every
      # time it is used.
      cat("\nDraws satisfying the signs: ", arw[["accepted"]], " of ",
          arw[["candidates"]], " (",
          format(round(100 * arw[["accepted"]] / arw[["candidates"]], 1), nsmall = 1),
          "%)\n", sep = "")
      cat("Effective sample size: ", arw[["effective_sample_size"]], " (",
          format(round(100 * arw[["effective_sample_size"]] / arw[["accepted"]], 1),
                 nsmall = 1),
          "% of those)\n", sep = "")
    }
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
      temp[, "Signif."] <- temp[, x[["model"]][["ci"]][1]] > 0 | temp[, x[["model"]][["ci"]][2]] < 0
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
    
    temp <- temp[which(lower.tri(matrix(1:(k * k), k), diag = TRUE)), , drop = FALSE]
    
    if (k == 1) {
      cat("\nVariance:\n\n")
    } else {
      cat("\nVariance-covariance matrix:\n\n") 
    }
    temp <- as.data.frame(temp)
    temp[, "Signif."] <- temp[, x[["model"]][["ci"]][1]] > 0 | temp[, x[["model"]][["ci"]][2]] < 0
    temp[, "Signif."] <- ifelse(temp[, "Signif."], "*", "")
    names(temp)[length(names(temp))] <- ""
    print(temp, ...)
  }

  .print_chains_summary(x[["chains"]])

  cat("\n")
  invisible(x)
}
