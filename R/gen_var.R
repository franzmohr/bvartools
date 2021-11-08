#' Vector Autoregressive Model Input
#' 
#' \code{gen_var} produces the input for the estimation of a vector autoregressive (VAR) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer vector of the lag order (default is \code{p = 2}).
#' @param exogen an optional time-series object of external regressors.
#' @param s an optional integer vector of the lag order of the external regressors (default is \code{s = 2}).
#' @param deterministic a character specifying which deterministic terms should
#' be included. Available values are \code{"none"}, \code{"const"} (default) for an intercept,
#' \code{"trend"} for a linear trend, and \code{"both"} for an intercept with a linear trend.
#' @param seasonal logical. If \code{TRUE}, seasonal dummy variables are
#' generated as additional deterministic terms. The amount of dummies depends on the frequency of the
#' time-series object provided in \code{data}.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param sv logical indicating whether time varying error variances should be estimated by
#' employing a stochastic volatility algorithm.
#' @param fcst integer. Number of observations saved for forecasting evaluation.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 5000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details The function produces the data matrices for vector autoregressive (VAR)
#' models, which can also include unmodelled, non-deterministic variables:
#' \deqn{A_0 y_t = \sum_{i=1}^{p} A_i y_{t - i} +
#' \sum_{i=0}^{s} B_i x_{t - i} +
#' C D_t + u_t,}
#' where
#' \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} coefficient matrix of contemporaneous endogenous variables,
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{x_t} is an M-dimensional vector of exogenous regressors and
#' \eqn{B_i} its corresponding \eqn{K \times M} coefficient matrix.
#' \eqn{D_t} is an N-dimensional vector of deterministic terms and
#' \eqn{C} its corresponding \eqn{K \times N} coefficient matrix.
#' \eqn{p} is the lag order of endogenous variables, \eqn{s} is the lag
#' order of exogenous variables, and \eqn{u_t} is an error term.
#' 
#' If an integer vector is provided as argument \code{p} or \code{s}, the function will
#' produce a distinct model for all possible combinations of those specifications.
#' 
#' If \code{tvp} is \code{TRUE}, the respective coefficients
#' of the above model are assumed to be time varying. If \code{sv} is \code{TRUE},
#' the error covariance matrix is assumed to be time varying.
#' 
#' @return An object of class \code{'bvarmodel'}, which contains the following elements:
#' \item{data}{A list of data objects, which can be used for posterior simulation. Element
#' \code{Y} is a time-series object of dependent variables. Element \code{Z} is a time-series
#' object of the regressors and element \code{SUR} is the corresponding matrix of regressors
#' in SUR form.}
#' \item{model}{A list of model specifications.}
#' 
#' @examples
#' 
#' # Load data
#' data("e1") 
#' e1 <- diff(log(e1))
#' 
#' # Generate model data
#' data <- gen_var(e1, p = 0:2, deterministic = "const")
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
gen_var <- function(data, p = 2, exogen = NULL, s = NULL,
                    deterministic = "const", seasonal = FALSE,
                    structural = FALSE, tvp = FALSE, sv = FALSE,
                    fcst = NULL,
                    iterations = 50000, burnin = 5000) {
  
  # Check data ----
  if (!"ts" %in% class(data)) {
    stop("Argument 'data' must be an object of class 'ts'.")
  }
  
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be an object of class 'ts'.")
    }
  }
  
  if (seasonal & !deterministic %in% c("const", "both")) {
    stop("Argument 'deterministic' must be either 'const' or 'both' when using 'seasonal = TRUE'.")
  }
  
  # Endogenous variables ----
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
  p_max <- max(p)
  
  model <- NULL
  model[["type"]] <- "VAR"
  model$endogen <- list("variables" = data_name,
                        "lags" = 0)
  
  temp <- data
  temp_name <- data_name
  if (p_max >= 1) {
    for (i in 1:p_max) {
      temp <- cbind(temp, stats::lag(data, -i))
      if (nchar(p_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0(data_name, ".", i_temp))
    }
  }
  
  # Exogenous variables ---- 
  if (!is.null(exogen)) {
    use_exo <- TRUE
    if (is.null(dimnames(exogen))) {
      tsp_temp <- stats::tsp(exogen)
      exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
      stats::tsp(exogen) <- tsp_temp
      dimnames(exogen)[[2]] <- "x"
    }
    exo_name <- dimnames(exogen)[[2]]
    m <- length(exo_name)
    s_max <- max(s)
    
    temp <- cbind(temp, exogen)
    if (nchar(s_max) > 2) {
      i_temp <- rep(0, nchar(s_max))
    } else {
      i_temp <- rep(0, 2)
    }
    i_temp <- paste0(i_temp, collapse = "")
    temp_name <- c(temp_name, paste0(exo_name, ".l", i_temp))
    if (s_max > 0) {
      for (i in 1:s_max) {
        temp <- cbind(temp, stats::lag(exogen, -i))
        if (nchar(s_max) > 2) {
          i_temp <- paste0(c(rep(0, nchar(s_max) - nchar(i)), i), collapse = "")
        } else {
          i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
        }
        temp_name <- c(temp_name, paste0(exo_name, ".l", i_temp))
      } 
    }
    
    model$exogen <- list("variables" = exo_name,
                         "lags" = 0)
  } else {
    use_exo <- FALSE
    s <- 0
    s_max <- 0
    m <- 0
  }
  
  temp <- stats::na.omit(temp)
  tt <- nrow(temp)
  det_name <- NULL
  
  if (deterministic %in% c("const", "both")) {
    temp <- cbind(temp, 1)
    temp_name <- c(temp_name, "const")
    det_name <- c(det_name, "const")
  }
  
  if (deterministic %in% c("trend", "both")) {
    temp <- cbind(temp, 1:tt)
    temp_name <- c(temp_name, "trend")
    det_name <- c(det_name, "trend")
  }
  
  if (seasonal) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(stats::cycle(temp) == 1)[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      for (i in 1:(freq - 1)) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = tt))
        temp_name <- c(temp_name, paste("season.", i, sep = ""))
        det_name <- c(det_name, paste("season.", i, sep = ""))
      }
    }
  }
  
  use_det <- FALSE
  if (length(det_name) > 0) {
    model[["deterministic"]] <- det_name
    use_det <- TRUE
  }
  
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
  
  # Is equal across all models
  y <- stats::ts(as.matrix(temp[, 1:k]), class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- stats::tsp(temp)
  dimnames(y)[[2]] <- temp_name[1:k]
  
  fcst_y <- NULL
  if (!is.null(fcst)) {
    fcst_y <- stats::window(y, start = stats::time(y)[nrow(y) - fcst + 1])
    y <- stats::window(y, end = stats::time(y)[nrow(y) - fcst])
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
      pos <- NULL
      model_i <- model
      if (i >= 1) {
        pos <- c(pos, k + 1:(k * i))
        model_i[["endogen"]][["lags"]] <- i
      }  
      if (use_exo) {
        pos <- c(pos, k + k * p_max + 1:(m * (j + 1)))
        model_i[["exogen"]][["lags"]] <- j
      }
      if (use_det) {
        pos <- c(pos, k + k * p_max + m * (s_max + 1) + 1:length(det_name))
      }
      
      x <- NULL
      z <- NULL
      if (length(pos) > 0) {
        x <- stats::ts(as.matrix(temp[, pos]), class = c("mts", "ts", "matrix")) 
        stats::tsp(x) <- stats::tsp(temp)
        dimnames(x)[[2]] <- temp_name[pos]
        
        if (!is.null(fcst)) {
          x <- stats::window(x, end = stats::time(x)[nrow(x) - fcst])
        }
        
        z <- kronecker(x, diag(1, k))
      }
      
      if (!is.null(y_A0)) {
        z <- cbind(z, y_A0)
      }
      
      result_i <- list("data" = list("Y" = y,
                                     "Z" = x,
                                     "SUR" = z,
                                     "TEST" = fcst_y),
                       "model" = model_i)
      
      result <- c(result, list(result_i)) 
      
    }
  }
  
  if (length(result) == 1) {
    result <- result[[1]]
  }
  
  class(result) <- append("bvarmodel", class(result)) 
  return(result)
}