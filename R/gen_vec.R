#' Vector Error Correction Model Input
#' 
#' \code{gen_vec} produces the input for the estimation of a vector error correction (VEC) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer vector of the lag order of the series in the (levels) VAR.
#' @param r an integer vector of the cointegration rank.
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer vector of the lag order of the exogenous variables of the series
#' in the (levels) VAR.
#' @param const a character specifying whether a constant term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param trend a character specifying whether a trend term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param seasonal a character specifying whether seasonal dummies should be included in the error
#' correction term (\code{"restricted"}) or in the non-cointegreation term as \code{"unrestricted"}
#' variables. If \code{NULL} (default) no seasonal terms will be added. The amount of dummy variables depends
#' on the frequency of the time-series object provided in \code{data}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param sv logical indicating whether time varying error variances should be estimated using
#' stochastic volatility algorithms.
# @param tvp a character vector specifying which parts of the model should be time varying.
# Possible elements are \code{"none"} (default), \code{"A0"}, \code{"alpha"}, \code{"beta"},
# \code{"Gamma"}, \code{"Upsilon"}, \code{"C"}, \code{"Sigma"}, \code{"all"}.
#' @param fcst integer. Number of observations saved for forecasting evaluation.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 5000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details The function produces the variable matrices of vector error correction (VEC)
#' models, which can also include exogenous variables:
#' \deqn{\Delta y_t = \Pi w_t + \sum_{i=1}^{p-1} \Gamma_{i} \Delta y_{t - i} + 
#' \sum_{i=0}^{s-1} \Upsilon_{i} \Delta x_{t - i} +
#' C^{UR} d^{UR}_t + u_t,}
#' where
#' \eqn{\Delta y_t} is a \eqn{K \times 1} vector of differenced endogenous variables,
#' \eqn{w_t} is a \eqn{(K + M + N^{R}) \times 1} vector of cointegration variables,
#' \eqn{\Pi} is a \eqn{K \times (K + M + N^{R})} matrix of cointegration parameters,
#' \eqn{\Gamma_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{\Delta x_t} is a \eqn{M \times 1} vector of differenced exogenous regressors,
#' \eqn{\Upsilon_i} is a \eqn{K \times M} coefficient matrix of exogenous regressors,
#' \eqn{d^{UR}_t} is a \eqn{N \times 1} vector of deterministic terms, and
#' \eqn{C^{UR}} is a \eqn{K \times N^{UR}} coefficient matrix of deterministic terms
#' that do not enter the cointegration term.
#' \eqn{p} is the lag order of endogenous variables and \eqn{s} is the lag
#' order of exogenous variables of the corresponding VAR model.
#' \eqn{u_t} is a \eqn{K \times 1} error term.
#' 
#' If an integer vector is provided as argument \code{p}, \code{s} or \code{r}, the function will
#' produce a distinct model for all possible combinations of those specifications.
#' 
#' If \code{tvp} is \code{TRUE}, the respective coefficients
#' of the above model are assumed to be time varying. If \code{sv} is \code{TRUE},
#' the diagonal elements of the error covariance matrix are assumed to be time varying.
#' 
#' @return An object of class \code{'bvecmodel'}, which contains the following elements:
#' \item{data}{A list of data objects, which can be used for posterior simulation. Element
#' \code{Y} is a time-series object of dependent variables. Element \code{W} is a timer-series
#' object of variables in the cointegration term and element \code{X} is a time-series
#' object of variables that do not enter the cointegration term. Element \code{SUR} contains a
#' matrix of element \code{X} in its SUR form.}
#' \item{model}{A list of model specifications.}
#' 
#' @examples 
#' 
#' # Load data
#' data("e6")
#' 
#' # Generate model data
#' data <- gen_vec(e6, p = 4, const = "unrestricted", season = "unrestricted")
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
gen_vec <- function(data, p = 2, exogen = NULL, s = 2, r = NULL,
                    const = NULL, trend = NULL, seasonal = NULL,
                    structural = FALSE, tvp = FALSE, sv = FALSE,
                    fcst = NULL,
                    iterations = 50000, burnin = 5000) {
  
  # rm(list = ls()[which(ls() != "data")]); p = 1:3; r = NULL; exogen = NULL; s = 1:2; const = "unrestricted"; trend = NULL; seasonal = "unrestricted"; structural = FALSE; iterations = 50000; burnin = 5000
  
  # Check data ----
  if (!"ts" %in% class(data)) {
    stop("Argument 'data' must be an object of class 'ts'.")
  }
  
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be an object of class 'ts'.")
    }
  }
  
  if (!is.null(const)) {
    if (!const %in% c("restricted", "unrestricted")) {
      stop("Specification for argument 'const' is not valid.")
    }
  }
  if (!is.null(trend)) {
    if (!trend %in% c("restricted", "unrestricted")) {
      stop("Specification for argument 'trend' is not valid.")
    }
  }
  if (!is.null(seasonal)) {
    if (!seasonal %in% c("restricted", "unrestricted")) {
      stop("Specification for argument 'seasonal' is not valid.")
    } else {
      if (is.null(const)) {
        stop("If argument 'seasonal' is specified, argument 'const' must be specified as well.")
      }
    }
  }
  
  if (0 %in% p) {
    stop("Argument 'p' must be at least 1 for any error correction model.")
  }
  if (0 %in% s) {
    stop("Argument 's' must be at least 1 for any error correction model.")
  }
  
  if (is.null(dimnames(data))) {
    tsp_temp <- stats::tsp(data)
    data <- stats::ts(as.matrix(data), class = c("mts", "ts", "matrix"))
    stats::tsp(data) <- tsp_temp
    dimnames(data)[[2]] <- "y"
  }
  
  if (NCOL(data) == 1 & structural) {
    stop("Model must contain at least two endogenous variables for structural analysis.")
  }
  
  data_name <- dimnames(data)[[2]]
  k <- NCOL(data)
  n_ect <- k
  p_max <- max(p)
  
  model <- NULL
  model[["type"]] <- "VEC"
  model[["endogen"]] <- list("variables" = dimnames(data)[[2]],
                             "lags" = 1)
  
  # Differenced endogenous variables
  diff_y <- diff(data)
  temp_name <- paste("d.", data_name, sep = "")
  temp <- diff_y
  
  # Endogenous ECT variables
  temp <- cbind(temp, stats::lag(data, -1))
  temp_name <- c(temp_name, paste("l.", data_name, sep = ""))
  
  # Exogenous ECT variables
  if (!is.null(exogen)) {
    use_exo <- TRUE
    if (is.null(dimnames(exogen))) {
      tsp_temp <- stats::tsp(exogen)
      exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
      stats::tsp(exogen) <- tsp_temp
      dimnames(exogen)[[2]] <- "x"
    }
    exog_name <- dimnames(exogen)[[2]]
    m <- length(exog_name)
    s_max <- max(s)
    
    # Non-deterministic exogenous ECT variables
    temp <- cbind(temp, stats::lag(exogen, -1))
    temp_name <- c(temp_name, paste("l.", exog_name, sep = ""))
    n_ect <- n_ect + m
    
    model[["exogen"]] <- list("variables" = dimnames(exogen)[[2]],
                              "lags" = 1)
  } else {
    use_exo <- FALSE
    s <- 0
    s_max <- 0
    m <- 0
  }
  
  # Lags of differenced endogenous variables
  if (p_max > 1) {
    for (i in 1:(p_max - 1)) {
      temp <- cbind(temp, stats::lag(diff_y, -i))
      if (nchar(p_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0("d.", data_name, ".l", i_temp))
    }
  }
  
  # Lags of differenecd exogenous variables
  if (use_exo) {
    # Add exogen s = 0
    diff_exog <- diff(exogen)
    temp <- cbind(temp, diff_exog)
    if (nchar(s_max) > 2) {
      i_temp <- rep(0, nchar(s_max))
    } else {
      i_temp <- rep(0, 2)
    }
    i_temp <- paste0(i_temp, collapse = "")
    temp_name <- c(temp_name, paste0("d.", exog_name, ".l", i_temp))
    if (s_max > 1) {
      for (i in 1:(s_max - 1)) {
        temp <- cbind(temp, stats::lag(diff_exog, -i))
        if (nchar(s_max) > 2) {
          i_temp <- paste0(c(rep(0, nchar(s_max) - nchar(i)), i), collapse = "")
        } else {
          i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
        }
        temp_name <- c(temp_name, paste0("d.", exog_name, ".l", i_temp))
      } 
    }
  }
  
  temp <- stats::na.omit(temp)
  ts_info <- stats::tsp(temp)
  
  # Final endogenous variables
  y <- stats::ts(as.matrix(temp[, 1:k]), class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- ts_info
  dimnames(y)[[2]] <- temp_name[1:k]
  
  # Forecast samples
  fcst_y <- NULL
  if (!is.null(fcst)) {
    fcst_y <- stats::window(y, start = stats::time(y)[nrow(y) - fcst + 1])
    y <- stats::window(y, end = stats::time(y)[nrow(y) - fcst])
    temp <- stats::window(temp, end = stats::time(temp)[nrow(temp) - fcst])
    ts_info <- stats::tsp(temp)
  }
  
  tt <- nrow(temp)
  
  ect <- matrix(temp[, k + 1:n_ect], tt)
  ect_names <- temp_name[k + 1:n_ect]
  
  x <- matrix(temp[, -(1:(k + n_ect))], tt)
  x_names <- temp_name[-(1:(k + n_ect))]
  
  det_name_r <- NULL
  det_name_ur <- NULL
  n_det_ur <- 0
  
  if (!is.null(const)) {
    if (const == "restricted") {
      ect <- cbind(ect, 1)
      ect_names <- c(ect_names, "const") 
      det_name_r <- c(det_name_r, "const") 
      n_ect <- n_ect + 1
    }
    
    if (const == "unrestricted") {
      x <- cbind(x, 1)
      x_names <- c(x_names, "const")
      det_name_ur <- c(det_name_ur, "const") 
      n_det_ur <- n_det_ur + 1
    }
  }
  
  if (!is.null(trend)) {
    if (trend == "restricted") {
      ect <- cbind(ect, 1:tt)
      ect_names <- c(ect_names, "trend")
      det_name_r <- c(det_name_r, "trend") 
      n_ect <- n_ect + 1
    }
    
    if (trend == "unrestricted") {
      x <- cbind(x, 1:tt)
      x_names <- c(x_names, "trend")
      det_name_ur <- c(det_name_ur, "trend")
      n_det_ur <- n_det_ur + 1
    }
  }
  
  if(!is.null(seasonal)) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(stats::cycle(temp) == 1)[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      seas <- NULL
      s_name <- NULL
      for (i in 1:(freq - 1)) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        seas <- cbind(seas, rep(s_temp, length.out = tt))
        s_name <- c(s_name, paste("season.", i, sep = ""))
      }
    }
    
    if (seasonal == "restricted") {
      ect <- cbind(ect, seas)
      ect_names <- c(ect_names, s_name)
      det_name_r <- c(det_name_r, s_name)
      n_ect <- n_ect + freq - 1
    }
    
    if (seasonal == "unrestricted") {
      x <- cbind(x, seas)
      x_names <- c(x_names, s_name)
      det_name_ur <- c(det_name_ur, s_name) 
      n_det_ur <- n_det_ur + length(s_name)
    }
  }
  
  use_det_r <- FALSE
  if (length(det_name_r) > 0) {
    use_det_r <- TRUE
    model[["deterministic"]][["restricted"]] <- det_name_r
  }
  use_det_ur <- FALSE
  if (length(det_name_ur) > 0) {
    use_det_ur <- TRUE
    model[["deterministic"]][["unrestricted"]] <- det_name_ur
  }
  
  if (is.null(r)) {
    if (n_ect > k) {
      r <- 0:k 
    } else {
      r <- 0:(k - 1)
    }
  } else {
    if (any(r > k)) {
      stop("Argument 'rank' must be smaller than or equal to the number of endogenous variables.")
    }
  }
  model[["rank"]] = 0
  
  if (class(structural) == "logical") {
    model[["structural"]] <- structural
  } else {
    stop("Argument 'structural' must be of class 'logical'.")
  }
  
  # TVP specifications ----
  if (class(tvp) == "logical") {
    model[["tvp"]] <- tvp 
  } else {
    stop("Argument 'tvp' must be of class 'logical'.")
  }
  if (class(sv) == "logical") {
    model[["sv"]] <- sv
  } else {
    stop("Argument 'sv' must be of class 'logical'.")
  }
  
  model[["iterations"]] <- iterations
  model[["burnin"]] <- burnin
  
  ect <- stats::ts(as.matrix(ect), class = c("mts", "ts", "matrix"))
  stats::tsp(ect) <- ts_info
  dimnames(ect)[[2]] <- ect_names
  
  if (length(x_names) > 0) {
    x <- stats::ts(as.matrix(x), class = c("mts", "ts", "matrix"))
    stats::tsp(x) <- ts_info
    dimnames(x)[[2]] <- x_names
  } else {
    x <- NULL
  }
  
  y_A0 <- NULL
  if (structural & k > 1) {
    y_A0 <- kronecker(-y, diag(1, k))
    pos <- NULL
    for (j in 1:k) {
      pos <- c(pos, (j - 1) * k + 1:j)
    }
    y_A0 <- y_A0[, -pos]
  }
  
  result <- NULL
  for (i in p) {
    for (j in s) {
      for (rank in r) {
        pos <- NULL
        model_i <- model
        if (i > 1) {
          pos <- c(pos, 1:(k * (i - 1)))
          model_i[["endogen"]][["lags"]] <- i
        }
        if (use_exo) {
          pos <- c(pos, k * (p_max - 1) + 1:(m * j))
          model_i[["exogen"]][["lags"]] <- j
        }
        if (use_det_ur) {
          pos <- c(pos, k * (p_max - 1) + m * s_max + 1:n_det_ur)
        }
        
        if (rank == 0 & length(pos) == 0 & !structural) {
          warning("Model with zero cointegration rank and no non-cointegration regressors is skipped.")
          next
        }
        
        model_i[["rank"]] = rank
        X <- NULL
        z <- NULL
        if (length(pos) > 0) {
          X <- stats::ts(as.matrix(x[, pos]), class = c("mts", "ts", "matrix")) 
          stats::tsp(X) <- ts_info
          dimnames(X)[[2]] <- x_names[pos]
        }
        
        z <- kronecker(cbind(ect, X), diag(1, k))
        
        if (!is.null(y_A0)) {
          z <- cbind(z, y_A0)
        }
        dimnames(z) <- NULL
        
        result_i <- list("data" = list("Y" = y,
                                       "W" = ect,
                                       "X" = X,
                                       "SUR" = z, 
                                       "TEST" = fcst_y),
                         "model" = model_i)
        
        result <- c(result, list(result_i)) 
      }
    }
  }
  
  if (!is.null(result)) {
    if (length(result) == 1) {
      result <- result[[1]]
    }
    class(result) <- append("bvecmodel", class(result))  
  } else {
    warning("Specification results in no output.")
  }
  return(result)
}