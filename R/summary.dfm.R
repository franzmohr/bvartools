#' Summarising Bayesian Dynamic Factor Models
#'
#' summary method for class \code{"dfm"}.
#'
#' @param object an object of class \code{"dfm"}, usually, a result of a call to
#' \code{\link{dfm}}.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.dfm} returns a list of class \code{"summary.dfm"},
#' which contains the following components:
#' \item{lambda}{A list of various summary statistics of the posterior
#' draws of the factor loadings.}
#' \item{factor}{A list of various summary statistics of the posterior
#' draws of the factors.}
#' \item{sigma_u}{A list of various summary statistics of the posterior
#' draws of the variance matrix of the measurement equation.}
#' \item{a}{A list of various summary statistics of the posterior
#' draws of the factor loadings.}
#' \item{sigma_v}{A list of various summary statistics of the posterior
#' draws of the variance matrix of the transition equation.}
#' \item{specifications}{a list containing information on the model specification.}
#'
#' @export
summary.dfm <- function(object, ci = .95, ...){

  tt <- NROW(object[["x"]])
  m <- NCOL(object[["x"]])
  n <- NCOL(object[["factor"]]) / tt
  n_lambda <- m * n
  p <- object$specifications$lags["p"]
  measure_names <- dimnames(object$x)[[2]]
  fac_names <- paste0("Factor.", 1:n)
  dim_names <- list(measure_names, fac_names)
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  
  # Measurement equation ----
  
  # Factor loadings
  temp <- summary(object[["lambda"]], quantiles = c(ci_low, .5, ci_high))
  means <- matrix(temp$statistics[, "Mean"], m)
  sds <- matrix(temp$statistics[, "SD"], m)
  naive_sd <- matrix(temp$statistics[, "Naive SE"], m)
  ts_sd <- matrix(temp$statistics[, "Time-series SE"], m)
  q_low <- matrix(temp$quantiles[, 1], m)
  median <- matrix(temp$quantiles[, 2], m)
  q_high <- matrix(temp$quantiles[, 3], m)

  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
  }
  
  result <- list(lambda = list(means = means,
                               median = median,
                               sd = sds,
                               naivesd = naive_sd,
                               tssd = ts_sd,
                               q_lower = q_low,
                               q_upper = q_high))
  
  # Factors
  means <- NULL
  median <- NULL
  sds <- NULL
  naive_sd <- NULL
  ts_sd <- NULL
  q_low <- NULL
  q_high <- NULL
  
  temp <- summary(object[["factor"]], quantiles = c(ci_low, .5, ci_high))
  for (i in 1:n) {
    if ("numeric" %in% class(temp$statistics)) {
      means <- cbind(means, matrix(temp$statistics["Mean"], tt))
      sds <- cbind(sds, matrix(temp$statistics["SD"], tt))
      naive_sd <- cbind(naive_sd, matrix(temp$statistics["Naive SE"], tt))
      ts_sd <- cbind(ts_sd, matrix(temp$statistics["Time-series SE"], tt))
      q_low <- cbind(q_low, matrix(temp$quantiles[1], tt))
      median <- cbind(median, matrix(temp$quantiles[2], tt))
      q_high <- cbind(q_high, matrix(temp$quantiles[3], tt)) 
    } else {
      means <- cbind(means, matrix(temp$statistics[i + n * 0:(tt - 1), "Mean"], tt))
      sds <- cbind(sds, matrix(temp$statistics[i + n * 0:(tt - 1), "SD"], tt))
      naive_sd <- cbind(naive_sd, matrix(temp$statistics[i + n * 0:(tt - 1), "Naive SE"], tt))
      ts_sd <- cbind(ts_sd, matrix(temp$statistics[i + n * 0:(tt - 1), "Time-series SE"], tt))
      q_low <- cbind(q_low, matrix(temp$quantiles[i + n * 0:(tt - 1), 1], tt))
      median <- cbind(median, matrix(temp$quantiles[i + n * 0:(tt - 1), 2], tt))
      q_high <- cbind(q_high, matrix(temp$quantiles[i + n * 0:(tt - 1), 3], tt)) 
    }
  }
  
  dim_names <- list(1:tt, fac_names)
  
  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
  }
  
  result[["factor"]] = list(means = means,
                            median = median,
                            sd = sds,
                            naivesd = naive_sd,
                            tssd = ts_sd,
                            q_lower = q_low,
                            q_upper = q_high)
  
  # Sigma u
  if (!is.null(object$sigma_u)) {
    temp <- summary(object$sigma_u, quantiles = c(ci_low, .5, ci_high))
    means <- matrix(temp$statistics[, "Mean"], m)
    sds <- matrix(temp$statistics[, "SD"], m)
    naive_sd <- matrix(temp$statistics[, "Naive SE"], m)
    ts_sd <- matrix(temp$statistics[, "Time-series SE"], m)
    q_low <- matrix(temp$quantiles[, 1], m)
    median <- matrix(temp$quantiles[, 2], m)
    q_high <- matrix(temp$quantiles[, 3], m)
  }
  
  dim_names <- list(measure_names, NULL)
  
  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
  }
  
  result[["sigma_u"]] = list(means = means,
                             median = median,
                             sd = sds,
                             naivesd = naive_sd,
                             tssd = ts_sd,
                             q_lower = q_low,
                             q_upper = q_high)

  # Transition equation ----
  # A
  if (p > 0) {
    temp <- summary(object[["a"]], quantiles = c(ci_low, .5, ci_high))
    if ("numeric" %in% class(temp$statistics)) {
      means <- matrix(temp$statistics["Mean"], n)
      sds <-matrix(temp$statistics["SD"], n)
      naive_sd <- matrix(temp$statistics["Naive SE"], n)
      ts_sd <- matrix(temp$statistics["Time-series SE"], n)
      q_low <- matrix(temp$quantiles[1], n)
      median <- matrix(temp$quantiles[2], n)
      q_high <- matrix(temp$quantiles[3], n)
    } else {
      means <- matrix(temp$statistics[, "Mean"], n)
      sds <- matrix(temp$statistics[, "SD"], n)
      naive_sd <- matrix(temp$statistics[, "Naive SE"], n)
      ts_sd <- matrix(temp$statistics[, "Time-series SE"], n)
      q_low <- matrix(temp$quantiles[, 1], n)
      median <- matrix(temp$quantiles[, 2], n)
      q_high <- matrix(temp$quantiles[, 3], n)
    } 
  }
  
  # Sigma v
  if (!is.null(object$sigma_v)) {
    temp <- summary(object$sigma_v, quantiles = c(ci_low, .5, ci_high))
    if (n == 1) {
      means <- matrix(temp$statistics["Mean"], n)
      sds <- matrix(temp$statistics["SD"], n)
      naive_sd <- matrix(temp$statistics["Naive SE"], n)
      ts_sd <- matrix(temp$statistics["Time-series SE"], n)
      q_low <- matrix(temp$quantiles[1], n)
      median <- matrix(temp$quantiles[2], n)
      q_high <- matrix(temp$quantiles[3], n)
    } else {
      means <- matrix(temp$statistics[, "Mean"], n)
      sds <- matrix(temp$statistics[, "SD"], n)
      naive_sd <- matrix(temp$statistics[, "Naive SE"], n)
      ts_sd <- matrix(temp$statistics[, "Time-series SE"], n)
      q_low <- matrix(temp$quantiles[, 1], n)
      median <- matrix(temp$quantiles[, 2], n)
      q_high <- matrix(temp$quantiles[, 3], n) 
    }
  }
  
  dim_names <- list(fac_names, NULL)
  
  if (!is.null(means)) {
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
  }
  
  result[["sigma_v"]] = list(means = means,
                             median = median,
                             sd = sds,
                             naivesd = naive_sd,
                             tssd = ts_sd,
                             q_lower = q_low,
                             q_upper = q_high)
  
  result$specifications <- object$specifications
  result$specifications$ci <- paste(c(ci_low, ci_high) * 100, "%", sep = "")
  
  class(result) <- "summary.dfm"
  return(result)
}
