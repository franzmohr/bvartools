#' Create Vector Error Correction Models
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
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param tvp logical indicating whether the model parameters are time varying.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param algorithm algorithm that should be used for posterior simulation. If \code{NULL}
#' (default), standard algorithms will be used. See 'Details' for available
#' non-standard options.
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
#' Available specifications for argument \code{varsel} are:
#' \itemize{
#'  \item{\code{"none"}: No variable selection algorithm is used.}
#'  \item{\code{"bvs"}: Bayesian variable selection as proposed in Korobilis (2013).}
#'  \item{\code{"ssvs"}: Stochastic search variable selection as proposed in George et al. (2008).}
#' }
#' 
#' Available specifications for argument \code{algorithm} are:
#' \itemize{
#'  \item{\code{"KLGS2010"}: Algorithm proposed in Koop, León-González & Strachan (2010).}
#' }
#' 
#' @return An object of class \code{'bvecmodel'}.
#' @examples 
#' 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_bvecmodel(e6, p = 4, r = 1,
#'                           const = "unrestricted", seasonal = "unrestricted",
#'                           iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' @references
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
create_bvecmodel <- function(data, p = 2, exogen = NULL, s = 2, r = NULL,
                             const = NULL, trend = NULL, seasonal = NULL,
                             structural = FALSE, error = "wishart", tvp = FALSE,
                             varsel = "none", algorithm = NULL,
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
  
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      stop("Invalid specification of argument 'error'.")
    }
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  # Exogenous variables ----
  use_exo <- !is.null(exogen)
  if (use_exo) {
    if (is.null(dimnames(exogen))) {
      # If 'exogen' is a simple ts object, transform it into a matrix object
      # to keep variable name information
      tsp_temp <- stats::tsp(exogen)
      exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
      stats::tsp(exogen) <- tsp_temp
      dimnames(exogen)[[2]] <- "exogen"
    }
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
  
  if (!varsel %in% c("none", "bvs", "ssvs")) {
    stop("Specification of argument 'varsel' is not supported.")
  }
  
  data_name <- dimnames(data)[[2]]
  k <- NCOL(data)
  n_ect <- k
  p_max <- max(p)
  
  if (is.null(algorithm)) {
    algo <- NULL
    if (tvp) {
      algo <- paste0(algo, "Tvp")
    } else {
      algo <- paste0(algo, "Normal")
    }
    if (error == "wishart") {
      algo <- paste0(algo, "Wishart")
    }
    if (error %in% c("gamma", "gamma+covar")) {
      algo <- paste0(algo, "Gamma")
    }
    if (error %in% c("sv", "sv+covar")) {
      algo <- paste0(algo, "Stochvol")
    }
    algo <- paste0("Vec", algo) 
  } else {
    if (!algorithm %in% c("KLGS2010")) {
      stop("Specified algorithm not recognized.")
    }
    if (algorithm == "KLGS2010") {
      algo <- "VecKlgs2010"
    }
  }
  
  model <- NULL
  model[["type"]] <- ifelse(length(data_name) == 1, "EC", "VEC")
  model[["algorithm"]] <- algo
  model[["k"]] <- k
  model[["p"]] <- 0L
  model[["m"]] <- 0L
  model[["s"]] <- 0L
  model[["n"]] <- 0L
  model[["n_restricted"]] <- 0L
  model[["rank"]] <- 0L
  model[["k_beta"]] <- 0L
  model[["varsel"]] <- varsel
  model[["endogen"]] <- dimnames(data)[[2]]
  if (use_exo) {
    exog_name <- dimnames(exogen)[[2]]
    model[["exogen"]] <- exog_name
  }
  
  # Differenced endogenous variables
  diff_y <- diff(data)
  temp_name <- paste("d.", data_name, sep = "")
  temp <- diff_y
  
  # Endogenous ECT variables
  temp <- cbind(temp, stats::lag(data, -1))
  temp_name <- c(temp_name, paste("l.", data_name, sep = ""))
  
  # Exogenous ECT variables
  if (!is.null(exogen)) {
    m <- length(exog_name)
    s_max <- max(s)
    
    # Non-deterministic exogenous ECT variables
    temp <- cbind(temp, stats::lag(exogen, -1))
    temp_name <- c(temp_name, paste("l.", exog_name, sep = ""))
    n_ect <- n_ect + m
    
    model[["m"]] <- m
  } else {
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
    model[["n"]] <- length(det_name_ur)
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
  
  if ("logical" %in% class(structural)) {
    model[["structural"]] <- structural
    if (structural) {
      model[["type"]] <- "SVECX" 
    }
  } else {
    stop("Argument 'structural' must be of class 'logical'.")
  }
  
  ## error ----
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      stop("Invalid specification of argument 'error'.")
    }
    model[["error"]] <- error
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  # TVP specifications ----
  if ("logical" %in% class(tvp)) {
    model[["tvp"]] <- tvp 
  } else {
    stop("Argument 'tvp' must be of class 'logical'.")
  }
  
  model[["iterations"]] <- as.integer(iterations)
  model[["burnin"]] <- as.integer(burnin)
  
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
  
  # Structural data
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
        
        # Build the matrix of regressors from object x
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
        
        model_i[["k_beta"]] <- as.integer(n_ect)
        model_i[["rank"]] = as.integer(rank)
        x_i <- NULL
        z <- NULL
        if (length(pos) > 0) {
          x_i <- stats::ts(as.matrix(x[, pos]), class = c("mts", "ts", "matrix")) 
          stats::tsp(x_i) <- ts_info
          dimnames(x_i)[[2]] <- x_names[pos]
          z <- kronecker(x_i, diag(1, k))
        }
        
        if (rank > 0) {
          z <- cbind(matrix(NA, tt * k, rank * k) , z)
        }
        
        if (!is.null(y_A0)) {
          z <- cbind(z, y_A0)
        }
        dimnames(z) <- NULL
        
        result_i <- list("model" = model_i,
                         "data" = list("original" = list("endogen" = data,
                                                         "exogen" = exogen),
                                       "train" = list("y" = y,
                                                      "w" = ect,
                                                      "x" = x_i,
                                                      "z" = z)))
        
        class(result_i) <- c("bvecmodel", "list") 
        
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