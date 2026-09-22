#' Prepare Forecast Input
#' 
#' Generates data matrices serving as input for forecasting simulation for objects
#' of class 'bvarmodel'.
#' 
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with deterministic data. If not
#' specified, the function will try to identify the deterministic terms
#' automatically. If this is not successful, an error message we be returned.
#' @param exogen a time-series object with the unmodelled, non-deterministic variables of the
#' model. Required if the model has such variables. See 'Details'.
#' @param ... additional arguments.
#'
#' @details The regressors of a forecast period contain the values of the unmodelled,
#' non-deterministic variables in that period and in the \code{s} periods before it. Argument
#' \code{exogen} therefore has to cover the last \code{s} periods of the estimation sample as well
#' as the \code{n_ahead} forecast periods, at the frequency of the model. A series that starts with
#' the first forecast period, or ends before the last one, is refused with a message that names the
#' periods it has to cover.
#'
#' If \code{deterministic} is not given, the deterministic terms are continued from the estimation
#' sample, which works for a constant, a linear trend and seasonal dummies.
#'
#' @return A list with elements \code{h}, the forecast horizon, and \code{x},
#' the out-of-sample regressors: \code{h} rows, one per period, by one column
#' per regressor. That is the compact layout, the same one a coefficient matrix
#' is \code{k} by; the SUR layout this used to return spread every regressor
#' over \code{k} columns and was \code{k^2} the size for no extra content.
#' The lags of the endogenous variables that the estimation sample does not
#' reach are the forecasts of earlier periods, which the forecast fills in as it
#' goes; they are zero here, since the samplers refuse a missing value anywhere
#' in \code{x}.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Generate forcast input
#' fcst_input <- prepare_forecast_input(model, n_ahead = 4)
#' 
#' 
#' @export
prepare_forecast_input.bvarmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL, ...) {
  
  # Input checks
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }
  
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be of class 'ts'.")
    }
  }
  
  if (!is.null(deterministic)) {
    if (!"ts" %in% class(deterministic)) {
      stop("Argument 'deterministic' must be of class 'ts'.")
    }
  }
  
  # Model specs
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  y <- object[["data"]][["train"]][["y"]]
  z <- object[["data"]][["train"]][["z"]]
  tt <- nrow(y)
  y_tsp <- stats::tsp(y)
  y_end <- y_tsp[2]
  y_freq <- y_tsp[3]
  
  if (m > 0 & is.null(exogen)) {
    stop("If parameters of unmodeled, non-deterministic variables are estimated, argument 'exogen' must be specified.")
  }
  
  if (n > 0 & is.null(deterministic)) {
    
    # Deducing deterministic terms from model input
    if (p > 0 | m > 0) {
      det_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]][-c(1:(k * p + m * (s + 1)))]
    } else {
      det_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]] 
    }
    
    deterministic <- matrix(NA, n_ahead + 1, n)
    dimnames(deterministic) <- list(NULL, det_names)
    
    ## constant
    if ("const" %in% det_names) {
      if (all(object[["data"]][["train"]][["x"]][, "const"] == 1)) {
        deterministic[, "const"] <- rep(1, n_ahead + 1)
      }
    }
    
    ## trend
    if ("trend" %in% det_names) {
      increment <- object[["data"]][["train"]][["x"]][tt, "trend"] - object[["data"]][["train"]][["x"]][tt - 1, "trend"]
      deterministic[, "trend"] <- seq(from = object[["data"]][["train"]][["x"]][tt, "trend"], by = increment, length.out = n_ahead + 1)
    }
    
    ## seasonal
    if (any(grepl("season", det_names, fixed = TRUE))) {
      # Detect seasonal pattern
      pos_season <- det_names[which(grepl("season", det_names, fixed = TRUE))]
      last_obs <- object[["data"]][["train"]][["x"]][tt, pos_season]
      for (i in 1:tt) {
        if (identical(last_obs, object[["data"]][["train"]][["x"]][i, pos_season])) {
          deterministic[, pos_season] <- object[["data"]][["train"]][["x"]][i + 0:n_ahead, pos_season]
          break
        }
      }
    }
    
    if (any(is.na(deterministic))) {
      stop("Could not identify all deterministic terms. Please specify argument 'deterministic' instead.") 
    }
    
    deterministic <- stats::ts(deterministic, start = y_end, frequency = y_freq)
  }
  
  n_tot <- k * p + m * (s + 1) + n
  
  pred_start <- stats::time(stats::ts(rep(NA, 2), start = y_end, frequency = y_freq))[-1]
  
  x <- NULL
  if (n_tot > 0) {
    x <- stats::ts(matrix(NA, n_ahead, n_tot), start = pred_start, frequency = y_freq)
    x_time <- stats::time(x)
    
    if (p > 0) {
      # The lags of the first forecast period are the last p periods of the
      # estimation sample, which is also where the forecast starts in time. They
      # used to be taken from the original series the model was created from:
      # every window of use_expanding_window() but the last, and a model cut
      # short by window(), then forecast from the end of the whole series.
      temp_p <- stats::embed(as.matrix(y), p)
      temp_p <- temp_p[nrow(temp_p),]
      for (i in 1:p) {
        if (i <= n_ahead) {
          x[i, ((i - 1) * k + 1):(p * k)] <- temp_p[1:((p - i + 1) * k)]
        }
      }

      # The lags the sample does not reach are forecasts of earlier horizons,
      # which the recursion writes in as it goes; nothing reads the value here.
      # It has to be finite all the same: the core refuses a NaN anywhere in
      # the forecast regressors, overwritten cells included, and so does
      # BayesTS for the model file write_to_hdf5() makes of this.
      lags <- x[, 1:(k * p), drop = FALSE]
      lags[is.na(lags)] <- 0
      x[, 1:(k * p)] <- lags
    }
    
    if (m > 0) {
      tsp_x <- stats::tsp(exogen)

      # The regressors of a forecast period are the exogenous values of that
      # period and the s before it, so the series has to reach back s periods
      # into the estimation sample and forward to the last forecast period.
      # stats::window() only warns about a start it cannot honour, and the
      # shorter matrix it returns then failed to fit into 'x' with R's "number
      # of items to replace", which said nothing about what was missing.
      need_start <- x_time[1] - s / y_freq
      need_end <- x_time[length(x_time)]
      tol <- 0.1 / y_freq
      if (abs(tsp_x[3] - y_freq) > 1e-8 || tsp_x[1] > need_start + tol || tsp_x[2] < need_end - tol) {
        stop("Argument 'exogen' must cover the periods from ", .ts_period_label(need_start, y_freq),
             " to ", .ts_period_label(need_end, y_freq), " at the frequency of the model (", y_freq,
             "): the ", s, " period(s) before the first forecast period, whose values enter its lagged ",
             "regressors, and the ", n_ahead, " forecast period(s). It covers ",
             .ts_period_label(tsp_x[1], tsp_x[3]), " to ", .ts_period_label(tsp_x[2], tsp_x[3]),
             " at frequency ", tsp_x[3], ".", call. = FALSE)
      }

      temp_x <- stats::ts(stats::embed(exogen, s + 1), end = tsp_x[2], frequency = tsp_x[3])
      temp_x <- stats::window(temp_x, start = x_time[1], end = x_time[length(x_time)])
      x[, k * p + 1:(m * (s + 1))] <- temp_x
    }
    
    if (n > 0) {
      deterministic <- stats::window(deterministic, start = x_time[1], end = x_time[length(x_time)])
      x[, k * p + m * (s + 1) + 1:n] <- deterministic
    }
    
    # Handed over as a plain matrix: the time series attributes describe the
    # periods the rows stand for and nothing downstream reads them, while the
    # sampler indexes rows by horizon.
    x <- matrix(as.numeric(x), n_ahead, n_tot)
  }
  
  result <- list("h" = as.integer(n_ahead),
                 "x" = x)
  
  return(result)
}

# Label of a period of a time series in the form in which it is passed to
# stats::window() or stats::ts(), e.g. c(2019, 4), which is how an error message
# asks for the periods a series has to cover.
.ts_period_label <- function(time, frequency) {
  year <- floor(time + 1e-6)
  if (frequency == 1) {
    return(as.character(year))
  }
  paste0("c(", year, ", ", round((time - year) * frequency) + 1, ")")
}
