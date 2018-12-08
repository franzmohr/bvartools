#' VAR Model Input
#' 
#' The function produces the input for the estimation of VAR model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p integer. The number of lags of endogenous variables. Defaults to 2.
#' @param exog an optional time-series object of exogenous variables.
#' @param q integer. The number of lags of the exogenous variables. Defaults to 2.
#' @param deterministic character specifying the used deterministic terms.
#' Possible values are "none", "const" (default), "trend" and "both".
#' @param seasonal logical specifying if seasonal dummies should be used. The amount of
#' dummies depends on the frequency of the time-series object `data`, i.e. `freqency(data) - 1`.
#' 
#' @return A list containing two elements:
#' \item{y}{A matrix of dependent variables.}
#' \item{x}{A matrix of regressor variables.}
#' 
#' @examples
#' data("e1")
#' data <- diff(log(e1))
#' 
#' var_input <- gen_var(data, p = 2)
#' 
#' @export
gen_var <- function(data, p = 2, exog = NULL, q = 2, deterministic = "const", seasonal = FALSE) {
  data_name <- dimnames(data)[[2]]
  temp_name <- data_name
  k <- NCOL(data)

  # Lags of endogenous variables
  temp <- data
  for (i in 1:p) {
    temp <- cbind(temp, stats::lag(data, -i))
    temp_name <- c(temp_name, paste(data_name, i, sep = "."))
  }
  
  # Lags of exogenous variables
  if (!is.null(exog)) {
    data_name <- dimnames(exog)[[2]]
    for (i in 0:q) {
      temp <- cbind(temp, stats::lag(exog, -i))
      temp_name <- c(temp_name, paste(data_name, i, sep = "."))
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
      for (i in 2:freq) {
        s_temp <- rep(0, freq)
        s_temp[i] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = t))
        temp_name <- c(temp_name, paste("season.", i - 1, sep = ""))
      }
    }
  }
  
  dimnames(temp)[[2]] <- temp_name
  
  y <- t(temp[, 1:k])
  x <- t(temp[, -(1:k)])
  result <- list("y" = y,
                 "x" = x)
}