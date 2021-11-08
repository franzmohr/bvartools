#' Summarising Bayesian VAR Coefficients
#'
#' summary method for class \code{"bvar"}.
#'
#' @param object an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param period integer. Index of the period, for which the summary statistics should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param x an object of class \code{"summary.bvar"}, usually, a result of a call to
#' \code{\link{summary.bvar}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.bvar} returns a list of class \code{"summary.bvar"},
#' which contains the following components:
#' \item{coefficients}{A list of various summary statistics of the posterior
#' draws of the VAR coefficients.}
#' \item{sigma}{A list of various summary statistics of the posterior
#' draws of the variance-covariance matrix.}
#' \item{specifications}{a list containing information on the model specification.}
#'
#' @export
summary.bvar <- function(object, ci = .95, period = NULL, ...){
  
  # Number of endogenous variables
  k <- object[["specifications"]][["dims"]][["K"]]
  tt <- NROW(object[["y"]])
  tvp <- object[["specifications"]][["tvp"]]
  if (any(unlist(tvp))) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
    period_long <- stats::time(object[["y"]])[period]
  } else {
    period_long <- NULL
  }
  
  # Obtain variable names
  dim_names <- list(NULL, NULL)
  if (!is.null(object[["y"]])) {
    dim_names[[1]] <- dimnames(object[["y"]])[[2]]
  } else {
    dim_names[[1]] <- paste("y", 1:k, sep = "")
  }
  
  # Extract names from data matrix
  m <- n <- o <- tot <- 0
  x_names <- NULL
  if (!is.null(object[["A"]])) {
    if (tvp[["A"]]) {
      m <- NCOL(object[["A"]][[1]]) / k
    } else {
      m <- NCOL(object[["A"]]) / k 
    }
    p <- m / k
    tot <- tot + m
    if (!is.null(object[["x"]])) {
      x_names <- c(x_names, dimnames(object[["x"]])[[2]][1:m])
    } else {
      for (i in 1:p) {
        x_names <- c(x_names, paste(dim_names[[1]], ".l", i, sep = ""))
      } 
    }
  }
  if (!is.null(object[["B"]])) {
    if (tvp[["B"]]) {
      n <- NCOL(object[["B"]][[1]]) / k
    } else {
      n <- NCOL(object[["B"]]) / k
    }
    tot <- tot + n
    if (!is.null(object[["x"]])) {
      x_names <- c(x_names, dimnames(object[["x"]])[[2]][m + 1:n])
    } else {
      x_names <- c(x_names, paste("x", 1:n, sep = ""))
    }
  }
  if (!is.null(object[["C"]])) {
    if (tvp[["C"]]) {
      o <- NCOL(object[["C"]][[1]]) / k
    } else {
      o <- NCOL(object[["C"]]) / k 
    }
    tot <- tot + o
    if (!is.null(object[["x"]])) {
      x_names <- c(x_names, dimnames(object[["x"]])[[2]][m + n + 1:o])
    } else {
      x_names <- c(x_names, paste("det", 1:m, sep = ""))
    }
  }
  if (!is.null(object[["A0"]])) {
    x_names <- c(x_names, dim_names[[1]])
  }
  dim_names[[2]] <- x_names
  
  # Non-error coefficients
  means <- NULL
  median <- NULL
  sds <- NULL
  naive_sd <- NULL
  ts_sd <- NULL
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  q_low <- NULL
  q_high <- NULL
  
  use_incl <- FALSE
  if (any(grepl("lambda", names(object)))) {
    use_incl <- TRUE
    incl <- NULL
  }
  
  vars <- c("A", "B", "C", "A0")
  for (i in vars) {
    if (!is.null(object[[i]])) {
      if (tvp[[i]]) {
        temp <- summary(object[[i]][[period]], quantiles = c(ci_low, .5, ci_high))
      } else {
        temp <- summary(object[[i]], quantiles = c(ci_low, .5, ci_high)) 
      }
      if ("numeric" %in% class(temp$statistics)) {
        means <- cbind(means, matrix(temp$statistics["Mean"], k))
        sds <- cbind(sds, matrix(temp$statistics["SD"], k))
        naive_sd <- cbind(naive_sd, matrix(temp$statistics["Naive SE"], k))
        ts_sd <- cbind(ts_sd, matrix(temp$statistics["Time-series SE"], k))
        q_low <- cbind(q_low, matrix(temp$quantiles[1], k))
        median <- cbind(median, matrix(temp$quantiles[2], k))
        q_high <- cbind(q_high, matrix(temp$quantiles[3], k)) 
      } else {
        means <- cbind(means, matrix(temp$statistics[, "Mean"], k))
        sds <- cbind(sds, matrix(temp$statistics[, "SD"], k))
        naive_sd <- cbind(naive_sd, matrix(temp$statistics[, "Naive SE"], k))
        ts_sd <- cbind(ts_sd, matrix(temp$statistics[, "Time-series SE"], k))
        q_low <- cbind(q_low, matrix(temp$quantiles[, 1], k))
        median <- cbind(median, matrix(temp$quantiles[, 2], k))
        q_high <- cbind(q_high, matrix(temp$quantiles[, 3], k)) 
      }
      if (use_incl) {
        var_temp <- paste0(i, "_lambda")
        if (var_temp %in% names(object)) {
          incl <- cbind(incl, matrix(colMeans(object[[var_temp]]), k))
        } else {
          incl <- cbind(incl, matrix(rep(NA_real_, ncol(object[[i]])), k))
        }
      }
    }
  }
  
  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
    if (use_incl) {
      dimnames(incl) <- dim_names
    }
  }
  
  result <- list(coefficients = list(means = means,
                                     median = median,
                                     sd = sds,
                                     naivesd = naive_sd,
                                     tssd = ts_sd,
                                     q_lower = q_low,
                                     q_upper = q_high))
  
  if (use_incl) {
    result[["coefficients"]][["lambda"]] = incl
  }
  
  dim_names <- list(dim_names[[1]], dim_names[[1]])
  
  # Error coefficients
  if (!is.null(object[["Sigma"]])) {
    if (tvp[["Sigma"]]) {
      temp <- summary(object[["Sigma"]][[period]], quantiles = c(ci_low, .5, ci_high))
    } else {
      temp <- summary(object[["Sigma"]], quantiles = c(ci_low, .5, ci_high)) 
    }
    if (k == 1) {
      means <- matrix(temp$statistics["Mean"], k)
      sds <- matrix(temp$statistics["SD"], k)
      naive_sd <- matrix(temp$statistics["Naive SE"], k)
      ts_sd <- matrix(temp$statistics["Time-series SE"], k)
      q_low <- matrix(temp$quantiles[1], k)
      median <- matrix(temp$quantiles[2], k)
      q_high <- matrix(temp$quantiles[3], k)
    } else {
      means <- matrix(temp$statistics[, "Mean"], k)
      sds <- matrix(temp$statistics[, "SD"], k)
      naive_sd <- matrix(temp$statistics[, "Naive SE"], k)
      ts_sd <- matrix(temp$statistics[, "Time-series SE"], k)
      q_low <- matrix(temp$quantiles[, 1], k)
      median <- matrix(temp$quantiles[, 2], k)
      q_high <- matrix(temp$quantiles[, 3], k)
    }

    
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
    
    result[["sigma"]] <- list(means = means,
                              median = median,
                              sd = sds,
                              naivesd = naive_sd,
                              tssd = ts_sd,
                              q_lower = q_low,
                              q_upper = q_high)
    
    if ("Sigma_lambda" %in% names(object)) {
      incl <- matrix(colMeans(object[["Sigma_lambda"]]), k)
      dimnames(incl) <- dim_names
      result[["sigma"]][["lambda"]] = incl
    }
  }
  
  result[["specifications"]] <- object[["specifications"]]
  result[["specifications"]][["ci"]] <- paste(c(ci_low, ci_high) * 100, "%", sep = "")
  result[["specifications"]][["period"]] <- period_long
  
  class(result) <- "summary.bvar"
  return(result)
}
