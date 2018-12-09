#' Vector Autoregressive Model Input
#' 
#' The function produces the input for the estimation of a vector autoregressive (VAR) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p integer for the lag order (default is `p = 2`).
#' @param exogen optionally, a time-series object of external regressors.
#' @param q optional integer for the lag order of the exogenous variables (default is `q = 2`).
#' @param deterministic type of deterministic regressors to include. Available values are
#' `"none"`, `"const"` (default), `"trend"`, and `"both"`.
#' @param seasonal logical. If `TRUE`, seasonal dummy variables will be generated. The
#' amount of dummies depends on the frequency of the time-series object provided in `data`.
#'  
#' @return A list containing two elements:
#' \item{y}{A matrix of the dependent variables.}
#' \item{x}{A matrix of the regressors.}
#' 
#' @examples
#' data("e1")
#' data <- diff(log(e1))
#' 
#' var_input <- gen_var(data, p = 2)
#' 
#' @export
gen_var <- function(data, p = 2, exogen = NULL, q = 2, deterministic = "const", seasonal = FALSE) {
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
    if (q >= 0) {
    for (i in 0:q) {
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
  result <- list("y" = y,
                 "x" = x)
}