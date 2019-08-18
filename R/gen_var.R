#' Vector Autoregressive Model Input
#' 
#' \code{gen_var} produces the input for the estimation of a vector autoregressive (VAR) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer of the lag order (default is \code{p = 2}).
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer of the lag order of the exogenous variables (default is \code{s = 2}).
#' @param deterministic a character specifying which deterministic terms should
#' be included. Available values are \code{"none"}, \code{"const"} (default),
#' \code{"trend"}, and \code{"both"}.
#' @param seasonal logical. If \code{TRUE}, seasonal dummy variables are
#' generated. The amount of dummies depends on the frequency of the
#' time-series object provided in \code{data}.
#' 
#' @details The function produces the variable matrices of a vector autoregressive (VAR)
#' model, which can also include exogenous variables:
#' \deqn{y_t = \sum_{i=1}^{p} A_i y_{t - i} +
#' \sum_{i=0}^{s} B_i x_{t - i} +
#' C D_t + u_t,}
#' where
#' \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{x_t} is an M-dimensional vector of exogenous regressors and
#' \eqn{B_i} its corresponding \eqn{K \times M} coefficient matrix.
#' \eqn{D_t} is an N-dimensional vector of deterministic terms and
#' \eqn{C} its corresponding \eqn{K \times N} coefficient matrix.
#' \eqn{p} is the lag order of endogenous variables, \eqn{s} is the lag
#' order of exogenous variables, and \eqn{u_t} is an error term.
#' 
#' In matrix notation the above model can be written as
#' \deqn{Y = P Z + U,}
#' where
#' \eqn{Y} is a \eqn{K \times T} matrix of endogenous variables,
#' \eqn{Z} is a \eqn{Kp + M(1+s) + N \times T} matrix of regressor variables,
#' and \eqn{U} is a \eqn{K \times T} matrix of errors. The function \code{gen_var}
#' generates the matrices \eqn{Y} and \eqn{Z}.
#' 
#' @return A list containing the following elements:
#' \item{Y}{a matrix of endogenous variables.}
#' \item{Z}{a matrix of regressor variables.}
#' 
#' @examples
#' data("e1")
#' e1 <- diff(log(e1))
#' data <- gen_var(e1, p = 2, deterministic = "const")
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
gen_var <- function(data, p = 2, exogen = NULL, s = 2, deterministic = "const", seasonal = FALSE) {
  if (!"ts" %in% class(data)) {
    stop("Data must be an object of class 'ts'.")
  }
  if (is.null(dimnames(data))) {
    tsp_temp <- stats::tsp(data)
    data <- stats::ts(as.matrix(data), class = c("mts", "ts", "matrix"))
    stats::tsp(data) <- tsp_temp
    dimnames(data)[[2]] <- "y"
  }
  data_name <- dimnames(data)[[2]]
  temp_name <- data_name
  k <- NCOL(data)
  
  model <- NULL
  model$endogenous <- list("variables" = dimnames(data)[[2]],
                           "lags" = 0)
  model$type <- "VAR"
  
  # Lags of endogenous variables
  temp <- data
  if (p >= 1) {
    for (i in 1:p) {
      temp <- cbind(temp, stats::lag(data, -i))
      temp_name <- c(temp_name, paste(data_name, i, sep = "."))
    }
    model$endogenous$lags <- p
  }
  
  # Lags of exogenous variables
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be an object of class 'ts'.")
    }
    if (is.null(dimnames(exogen))) {
      tsp_temp <- stats::tsp(exogen)
      exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
      stats::tsp(exogen) <- tsp_temp
      dimnames(exogen)[[2]] <- "x"
    }
    data_name <- dimnames(exogen)[[2]]
    if (s >= 0) {
      for (i in 0:s) {
        temp <- cbind(temp, stats::lag(exogen, -i))
        temp_name <- c(temp_name, paste(data_name, i, sep = "."))
      }
    }
    model$exogen <- list("variables" = dimnames(exogen)[[2]],
                         "lags" = s)
  }
  
  temp <- stats::na.omit(temp)
  t <- nrow(temp)
  det_name <- NULL
  
  if (deterministic %in% c("const", "both")) {
    temp <- cbind(temp, 1)
    temp_name <- c(temp_name, "const")
    det_name <- c(det_name, "const")
  }
  
  if (deterministic %in% c("trend", "both")) {
    temp <- cbind(temp, 1:t)
    temp_name <- c(temp_name, "trend")
    det_name <- c(det_name, "trend")
  }
  
  if (seasonal) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(floor(stats::time(temp)) == stats::time(temp))[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      for (i in 1:(freq - 1)) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = t))
        temp_name <- c(temp_name, paste("season.", i, sep = ""))
        det_name <- c(det_name, paste("season.", i, sep = ""))
      }
    }
  }
  
  if (length(det_name) > 0) {
    model$deterministic <- det_name
  }
  
  y <- matrix(t(temp[, 1:k]), k, dimnames = list(temp_name[1:k], NULL))
  attr(y, "ts_info") <- stats::tsp(temp)
  
  x <- matrix(t(temp[, -(1:k)]), length(temp_name) - k,
              dimnames = list(temp_name[-(1:k)], NULL))
  
  result <- list("Y" = y,
                 "Z" = x,
                 "model" = model)
  
  class(result) <- append("bvarmodel", class(result))
  return(result)
}