#' Create a Vector Error Correction Model
#' 
#' Produces the input for the estimation of a vector error correction (VEC) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer vector of the lag order of the series in the (levels) VAR. Thus, the
#' resulting model's lag will be \eqn{p - 1}. See 'Details'.
#' @param r an integer vector of the cointegration rank. See 'Details'.
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer vector of the lag order of the exogenous variables of the series
#' in the (levels) VAR. Thus, the resulting model's lag will be \eqn{s - 1}. See 'Details'.
#' @param const a character specifying whether a constant term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param trend a character specifying whether a trend term enters the error correction
#' term (\code{"restricted"}) or the non-cointegration term as an \code{"unrestricted"} variable.
#' If \code{NULL} (default) no constant term will be added.
#' @param seasonal a character specifying whether seasonal dummies should be included in the error
#' correction term (\code{"restricted"}) or in the non-cointegreation term as \code{"unrestricted"}
#' variables. If \code{NULL} (default) no seasonal terms will be added. The amount of dummy variables
#' will be automatically detected and depends on the frequency of the time-series object provided
#' in \code{data}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
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
#' \eqn{d^{UR}_t} is a \eqn{N^{UR} \times 1} vector of deterministic terms, and
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
#' of the above model are assumed to be time varying.
#' 
#' Argument \code{error} specifies the structure of the covariance matrix of
#' the error term and how it is estimated. Possible specifications are:
#' \itemize{
#'  \item{\code{"wishart"}: The covariance is estimated using a Wishart prior.}
#'  \item{\code{"gamma"}: Only the diagonal elements of the covariance matrix are estimated using a gamma prior.
#' Off-diagonal elements are not estimated and set to zero.}
#'  \item{\code{"gamma+covar"}: The diagonal elements of the covariance matrix are estimated using a gamma prior.
#' Covariances are estimated based on a triangular decomposition.}
#'  \item{\code{"sv"}: Only the diagonal elements of the covariance matrix are estimated using a stochastic volatility
#' algorithm. Off-diagonal elements are not estimated and set to zero.}
#'  \item{\code{"sv+covar"}: Only the diagonal elements of the covariance matrix are estimated using a stochastic volatility
#' algorithm. Covariances are estimated based on a triangular decomposition.}
#' }
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
#' # Create model
#' model <- create_vec_model(e6, p = 4, r = 1,
#'                           const = "unrestricted", seasonal = "unrestricted",
#'                           iterations = 100, burnin = 10)
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
create_vec_model <- function(data, p = 2, exogen = NULL, s = 2, r = NULL,
                             const = NULL, trend = NULL, seasonal = NULL,
                             structural = FALSE, error = "wishart", tvp = FALSE,
                             iterations = 20000, burnin = 2000) {
  
  # Input checks ----
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
  
  ## errors ----
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      stop("Invalid specification of argument 'error'.")
    }
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  if (is.null(dimnames(data))) {
    tsp_temp <- stats::tsp(data)
    data <- stats::ts(as.matrix(data), class = c("mts", "ts", "matrix"))
    stats::tsp(data) <- tsp_temp
    dimnames(data)[[2]] <- "y"
  }
  
  if (NCOL(data) == 1 & structural) {
    structural <- FALSE
    if (error == "gamma+covar") {
      error <- "gamma"
    }
    if (error == "sv+covar") {
      error <- "sv"
    }
  }
  
  if (structural & error %in% c("wishart", "gamma+covar", "sv+covar")) {
    stop(paste0("Structural models cannot be estimated with argument 'error' specified as '", error,"'."))
  }
  
  data_name <- dimnames(data)[[2]]
  k <- NCOL(data)
  n_ect <- k
  p_max <- max(p)
  
  model <- NULL
  model[["type"]] <- "VEC"
  model[["k"]] <- k
  model[["p"]] <- 0
  model[["m"]] <- 0
  model[["s"]] <- 0
  model[["n_restricted"]] <- 0
  model[["n_unrestricted"]] <- 0
  
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
    
    model$m <- m
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
    model[["n_restricted"]] <- length(det_name_r)
  }
  use_det_ur <- FALSE
  if (length(det_name_ur) > 0) {
    use_det_ur <- TRUE
    model[["n_unrestricted"]] <- length(det_name_ur)
  }
  
  if (is.null(r)) {
    if (n_ect > k) {
      r <- 0:k 
    } else {
      r <- 0:(k - 1)
    }
    message("Argument rank 'r' not specified. Generating models for r = ", paste0(r, collapse = ", "), ".")
  } else {
    if (any(r > k)) {
      stop("Argument rank 'r' must be smaller than or equal to the number of endogenous variables.")
    }
  }
  model[["rank"]] <- 0
  
  if ("logical" %in% class(structural)) {
    model[["structural"]] <- structural
  } else {
    stop("Argument 'structural' must be of class 'logical'.")
  }
  
  ## error ----
  model[["error"]] <- error
  
  # TVP specifications ----
  if ("logical" %in% class(tvp)) {
    model[["tvp"]] <- tvp 
  } else {
    stop("Argument 'tvp' must be of class 'logical'.")
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
        }
        if (i >= 1) {
          model_i[["p"]] <- i
        }
        if (use_exo) {
          pos <- c(pos, k * (p_max - 1) + 1:(m * j))
          model_i[["s"]] <- j
        }
        if (use_det_ur) {
          pos <- c(pos, k * (p_max - 1) + m * s_max + 1:n_det_ur)
        }
        
        # if (rank == 0 & length(pos) == 0 & !structural & !error %in% c("gamma+covar", "sv+covar")) {
        #   warning("Model with zero cointegration rank and no non-cointegration regressors is skipped.")
        #   next
        # }
        
        model_i[["rank"]] = rank
        X <- NULL
        z <- NULL
        if (length(pos) > 0) {
          X <- stats::ts(as.matrix(x[, pos]), class = c("mts", "ts", "matrix")) 
          stats::tsp(X) <- ts_info
          dimnames(X)[[2]] <- x_names[pos]
          z <- kronecker(X, diag(1, k))
        }
        
        if (rank > 0) {
          z <- cbind(matrix(NA, tt * k, rank * k) , z)
        }
        
        if (!is.null(y_A0)) {
          z <- cbind(z, y_A0)
        }
        dimnames(z) <- NULL
        
        result_i <- list("data" = list("y" = y,
                                       "w" = ect,
                                       "x" = X,
                                       "z" = z),
                         "model" = model_i)
        
        class(result_i) <- append("bvecmodel", class(result_i)) 
        
        result <- c(result, list(result_i)) 
      }
    }
  }
  
  if (!is.null(result)) {
    if (length(result) == 1) {
      result <- result[[1]]
    } else {
      class(result) <- append("modellist", class(result)) 
    }
    
  } else {
    warning("Specification results in no output.")
  }
  return(result)
}