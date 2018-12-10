#' Vector Error Correction Model Input
#' 
#' The function produces the input for the estimation of a vector error correction (VEC) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer for the lag order of the series (levels) in the VAR (default is `p = 2`).
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer for the lag order of the exogenous variables of the series
#' (levels) in the VAR (default is `q = 2`).
#' @param const a character specifying whether a constant term enters the error correction term
#' (`"restricted"`) or is `"unrestricted"`. If `NULL` (default) no constant term will be added.
#' @param trend a character specifying whether a linear trend term enters the error correction term
#' (`"restricted"`) or is `"unrestricted"`. If `NULL` (default) no trend term will be added.
#' @param seasonal a character specifying whether seasonal dummies should be included in the error
#' correction term (`"restricted"`) or as `"unrestricted"` variables. If `NULL` (default) no
#' seasonal terms will be added. The amount of dummies depends on the frequency of the time-series
#' object provided in `data`. 
#' 
#' @details The function produces the input matrices for the estimation of
#' the vector error correction (VEC) model
#' \deqn{\Delta y_t =  \Pi ECT_t + \sum_{i=1}^{p-1} \Gamma_i \Delta y_{t - i} + \sum_{i=0}^{s-1} \Upsilon_i \Delta x_{t - i} + C D_t + u_t,}
#' where
#' \eqn{\Delta y_t} is a \eqn{K \times 1} vector of differenced endogenous variables,
#' \eqn{ECT_t} is a \eqn{K + M + N_{co} \times 1} vector of cointegration variables,
#' \eqn{\Pi} is a \eqn{K \times K + M + N_{co}} matrix of cointegration parameters,
#' \eqn{\Gamma_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{\Delta x_t} is a \eqn{M \times 1} vector of differenced exogenous regressors,
#' \eqn{\Upsilon_i} is a \eqn{K \times M} coefficient matrix of exogenous regressors,
#' \eqn{d_t} is a \eqn{N \times 1} vector of deterministic terms, and
#' \eqn{D} is a \eqn{K \times N} coefficient matrix of deterministic terms.
#' \eqn{p} is the lag order of endogenous variables, \eqn{s} is the lag
#' order of exogenous variables, and \eqn{u_t} is a \eqn{K \times 1} error term.
#' 
#' In matrix notation the above model can be written as
#' \deqn{Y = \Pi ECT + \Gamma X + U,}
#' where
#' \eqn{Y} is a \eqn{K \times T} matrix of differenced endogenous variables,
#' \eqn{ECT} is a \eqn{K + M + N_{co} \times T} matrix of cointegration
#' variables, \eqn{X} is a \eqn{K(p - 1) + Ms + N \times T} matrix of
#' differenced regressor variables and unrestricted deterministic terms.
#' \eqn{U} is a \eqn{K \times T} matrix of errors. The function `gen_vec`
#' generates the matrices \eqn{Y}, \eqn{ECT} and \eqn{X}.
#' 
#' @return A list containing the following elements
#' \item{Y}{A matrix of differenced dependent variables.}
#' \item{ECT}{A matrix of regressor variables in levels and restricted
#' deterministic terms.}
#' \item{X}{A matrix of differenced regressor variables and unrestricted
#' deterministic terms.}
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' 
#' @examples
#' data("e6")
#' 
#' vec_input <- gen_vec(e6, p = 4, const = "unrestricted", seasonal = "unrestricted")
#' 
#' @export
gen_vec <- function(data, p = 2, exogen = NULL, s = 2, const = NULL, trend = NULL, seasonal = NULL) {
  data_name <- dimnames(data)[[2]]
  k <- NCOL(data)
  
  diff_y <- diff(data)
  temp_name <- paste("d.", data_name, sep = "")
  temp <- diff_y
  
  temp <- cbind(temp, stats::lag(data, -1))
  temp_name <- c(temp_name, paste("l.", data_name, sep = ""))
  n_ect <- k
  
  if (!is.null(exogen)) {
    exog_name <- dimnames(exogen)[[2]]
    temp <- cbind(temp, stats::lag(exogen, -1))
    temp_name <- c(temp_name, paste("l.", exog_name, sep = ""))
    n_ect <- n_ect + NCOL(exogen)
  }
  
  # Lags of differenced endogenous variables
  if (p > 1) {
    for (i in 1:(p - 1)) {
      temp <- cbind(temp, stats::lag(diff_y, -i))
      temp_name <- c(temp_name, paste("d", data_name, i, sep = "."))
    } 
  }
  
  # Lags of exogenous variables
  if (!is.null(exogen)) {
    diff_exog <- diff(exogen)
    temp <- cbind(temp, diff_exog)
    temp_name <- c(temp_name, paste("d", exog_name, 0, sep = "."))
    if (s >= 2) {
      for (i in 1:(s - 1)) {
        temp <- cbind(temp, stats::lag(diff_exog, -i))
        temp_name <- c(temp_name, paste("d", exog_name, i, sep = "."))
      } 
    }
  }
  
  temp <- stats::na.omit(temp)
  t <- nrow(temp)
  t_info <- stats::tsp(temp)
  
  y <- matrix(temp[, 1:k], t)
  y_names <- temp_name[1:k]
  ect <- matrix(temp[, k + 1:n_ect], t)
  ect_names <- temp_name[k + 1:n_ect]
  x <- matrix(temp[, -(1:(k + n_ect))], t)
  x_names <- temp_name[-(1:(k + n_ect))]
  
  if (!is.null(const)) {
    if (const == "restricted") {
      ect <- cbind(ect, 1)
      ect_names <- c(ect_names, "const") 
      n_ect <- n_ect + 1
    }
    
    if (const == "unrestricted") {
      x <- cbind(x, 1)
      x_names <- c(x_names, "const") 
    }
  }
  
  if (!is.null(trend)) {
    if (trend == "restricted") {
      ect <- cbind(ect, 1:t)
      ect_names <- c(ect_names, "trend")
      n_ect <- n_ect + 1
    }
    
    if (trend == "unrestricted") {
      x <- cbind(x, 1:t)
      x_names <- c(x_names, "trend")
    }
  }
  
  if(!is.null(seasonal)) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(floor(stats::time(temp)) == stats::time(temp))[1]
      pos <- rep(1:4, 2)[pos:(pos + (freq - 2))]
      seas <- NULL
      s_name <- NULL
      for (i in 1:3) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        seas <- cbind(seas, rep(s_temp, length.out = t))
        s_name <- c(s_name, paste("season.", i, sep = ""))
      }
    }
    
    if (seasonal == "restricted") {
      ect <- cbind(ect, seas)
      ect_names <- c(ect_names, s_name)
      n_ect <- n_ect + freq - 1
    }
    
    if (seasonal == "unrestricted") {
      x <- cbind(x, seas)
      x_names <- c(x_names, s_name)
    }
  }
  
  temp <- cbind(y, ect, x)
  dimnames(temp)[[2]] <- c(y_names, ect_names, x_names)
  
  y <- t(temp[, 1:k])
  ect <- t(temp[, k + 1:n_ect])
  x <- t(temp[, -(1:(k + n_ect))])
  result <- list("Y" = y,
                 "ECT" = ect,
                 "X" = x)
}