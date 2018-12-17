#' Vector Autoregressive Model Input
#' 
#' The function produces the input for the estimation of a vector autoregressive (VAR) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer for the lag order (default is `p = 2`).
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer for the lag order of the exogenous variables (default is `q = 2`).
#' @param deterministic a character specifying type of deterministic
#' regressors to include. Available values are `"none"`, `"const"` (default),
#' `"trend"`, and `"both"`.
#' @param seasonal logical. If `TRUE`, seasonal dummy variables will be
#' generated. The amount of dummies depends on the frequency of the
#' time-series object provided in `data`.
#' 
#' @details The function produces the input matrices for the estimation of the model
#' \deqn{y_t = \sum_{i=1}^{p} A_i y_{t - i} + \sum_{i=0}^{s} B_i x_{t - i} + C D_t + u_t,}
#' where
#' \eqn{y_t} is a \eqn{K \times 1} vector of endogenous variables,
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{x_t} is a \eqn{M \times 1} vector of exogenous regressors,
#' \eqn{B_i} is a \eqn{K \times M} coefficient matrix of exogenous regressors,
#' \eqn{D_t} is a \eqn{N \times 1} vector of deterministic terms, and
#' \eqn{C} is a \eqn{K \times N} coefficient matrix of deterministic terms.
#' \eqn{p} is the lag order of endogenous variables, \eqn{s} is the lag
#' order of exogenous variables, and \eqn{u_t} is a \eqn{K \times 1} error term.
#' 
#' In matrix notation the above model can be written as
#' \deqn{Y = A Z + U,}
#' where
#' \eqn{Y} is a \eqn{K \times T} matrix of endogenous variables,
#' \eqn{Z} is a \eqn{Kp \times M(1+s) + N x T} matrix of regressor variables,
#' and \eqn{U} is a \eqn{K \times T} matrix of errors. The function `gen_var`
#' generates the matrices \eqn{Y} and \eqn{Z}.
#' 
#' @return A list containing the following elements
#' \item{Y}{A matrix of dependent variables.}
#' \item{Z}{A matrix of regressor variables.}
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' @examples
#' data("e1")
#' data <- diff(log(e1))
#' 
#' var_input <- gen_var(data, p = 2)
#' 
#' @export
gen_var <- function(data, p = 2, exogen = NULL, s = 2, deterministic = "const", seasonal = FALSE) {
  data_name <- dimnames(data)[[2]]
  temp_name <- data_name
  k <- NCOL(data)

  # Lags of endogenous variables
  temp <- data
  if (p >= 1) {
    for (i in 1:p) {
      temp <- cbind(temp, stats::lag(data, -i))
      temp_name <- c(temp_name, paste(data_name, i, sep = "."))
    }
  }

  # Lags of exogenous variables
  if (!is.null(exogen)) {
    data_name <- dimnames(exogen)[[2]]
    if (s >= 0) {
    for (i in 0:s) {
      temp <- cbind(temp, stats::lag(exogen, -i))
      temp_name <- c(temp_name, paste(data_name, i, sep = "."))
    }
    }
  }
  
  temp <- stats::na.omit(temp)
  t <- nrow(temp)
  
  if (deterministic %in% c("const", "both")) {
    temp <- cbind(temp, 1)
    temp_name <- c(temp_name, "const") 
  }
  
  if (deterministic %in% c("trend", "both")) {
    temp <- cbind(temp, 1:t)
    temp_name <- c(temp_name, "trend") 
  }
  
  if (seasonal) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(floor(stats::time(temp)) == stats::time(temp))[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      for (i in 1:3) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = t))
        temp_name <- c(temp_name, paste("season.", i, sep = ""))
      }
    }
  }
  
  dimnames(temp)[[2]] <- temp_name
  
  y <- t(temp[, 1:k])
  x <- t(temp[, -(1:k)])
  result <- list("Y" = y,
                 "Z" = x)
  return(result)
}