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
#' @param y_start an optional time-series object with the values of the endogenous variables
#' the forecast starts from. The forecast periods are the \code{n_ahead} periods after its last
#' observation, which may lie after the end of the estimation sample. If \code{NULL}
#' (default), the forecast starts after the end of the estimation sample. See 'Details'.
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
#' sample, which works for a constant, a linear trend, seasonal dummies and the dummy variables of
#' \code{\link{add_dummy_variables}}.
#'
#' \strong{Starting a forecast after the end of the estimation sample.} By default the first
#' forecast period is the one after the last period of the estimation sample, and its lagged
#' regressors are the last \code{p} observations of that sample. Argument \code{y_start} moves
#' the start: the first forecast period is then the one after the last observation of
#' \code{y_start}, and the lags are its last \code{p} observations. This is how a forecast
#' starts from a value that was not part of the estimation, such as a nowcast of the current
#' month or a first release that came in after the model was estimated, without estimating the
#' model again. \code{y_start} holds the endogenous variables of the model, by name and in the
#' units of the estimation sample, at the frequency of the model. It has to end in the last
#' period of the estimation sample or after it, and its last \code{p} observations must not be
#' missing; values before those are not read. When it ends in the last period of the estimation
#' sample and holds the same values, the result is the same as without it.
#'
#' The other inputs follow the new forecast periods: \code{exogen} has to cover the \code{s}
#' periods before the first of them and the \code{n_ahead} periods themselves, and
#' \code{deterministic} the forecast periods. The automatic continuation of the deterministic
#' terms counts through the periods between the end of the sample and the start of the
#' forecast, so that a trend keeps counting and seasonal dummies fall on the calendar months or
#' quarters of the forecast periods.
#'
#' The coefficients, and for models with time varying coefficients or stochastic volatility
#' the states, are those estimated on the sample. Nothing is learnt from the start values
#' beyond their role as lags. Some uses of a forecast therefore need it to start right after
#' the sample: see \code{\link{add_posterior_forecasts}} and
#' \code{\link{add_predictive_loglik}}. A model estimated with \code{missing = "estimate"}
#' starts each draw's forecast from its own completed panel and does not take \code{y_start}.
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
#' If \code{y_start} is given, the list has two more elements: \code{start}, the time of
#' the first forecast period, and \code{y}, the observations of \code{y_start} after the end of
#' the estimation sample together with the last \code{p} the forecast starts from, a time-series
#' object with one column per endogenous variable.
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
#' # Start the forecast two quarters after the end of the sample, from
#' # values the model was not estimated with
#' short <- create_bvarmodel(window(e1, end = c(1978, 4)), p = 2,
#'                           deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' fcst_input <- prepare_forecast_input(short, n_ahead = 4,
#'                                      y_start = window(e1, end = c(1979, 2)))
#' fcst_input[["start"]]
#' 
#' 
#' @export
prepare_forecast_input.bvarmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL,
                                             y_start = NULL, ...) {
  
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
  
  # Where the forecast starts: 'gap' periods after the end of the sample, none
  # without start values, which leaves every line below as it was.
  gap <- 0L
  if (!is.null(y_start)) {
    start_values <- .check_forecast_start_values(object, y_start)
    y_start <- start_values[["y"]]
    gap <- start_values[["gap"]]
  }
  # Periods from the end of the sample to the last forecast period
  n_det <- gap + n_ahead

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
    
    # One row for the last period of the sample and one for every period after
    # it up to the last forecast period. With start values after the end of the
    # sample that includes the periods in between, which a trend and seasonal
    # dummies have to be counted through.
    deterministic <- matrix(NA, n_det + 1, n)
    dimnames(deterministic) <- list(NULL, det_names)

    ## constant
    if ("const" %in% det_names) {
      if (all(object[["data"]][["train"]][["x"]][, "const"] == 1)) {
        deterministic[, "const"] <- rep(1, n_det + 1)
      }
    }

    ## trend
    if ("trend" %in% det_names) {
      increment <- object[["data"]][["train"]][["x"]][tt, "trend"] - object[["data"]][["train"]][["x"]][tt - 1, "trend"]
      deterministic[, "trend"] <- seq(from = object[["data"]][["train"]][["x"]][tt, "trend"], by = increment, length.out = n_det + 1)
    }

    ## seasonal
    if (any(grepl("season", det_names, fixed = TRUE))) {
      pos_season <- det_names[which(grepl("season", det_names, fixed = TRUE))]
      deterministic[, pos_season] <-
        .continue_seasonal_pattern(object[["data"]][["train"]][["x"]][, pos_season, drop = FALSE], n_det)
    }

    ## dummy variables of add_dummy_variables(), continued by the rule each
    ## was built by rather than by pattern. Only in the periods the forecast
    ## reads, so that a dummy given as data need not cover the periods between
    ## the end of the sample and the start of the forecast.
    spec <- object[["model"]][["dummy_variables"]]
    if (!is.null(spec)) {
      det_times <- y_end + (gap + 0:n_ahead) / y_freq
      deterministic[gap + 0:n_ahead + 1, spec[["name"]]] <- .continue_dummy_variables(object, det_times)
    }

    # The periods the forecast reads: the one before its first period and the
    # forecast periods, as without start values
    deterministic <- deterministic[gap + 0:n_ahead + 1, , drop = FALSE]

    if (any(is.na(deterministic))) {
      stop("Could not identify all deterministic terms. Please specify argument 'deterministic' instead.")
    }

    deterministic <- stats::ts(deterministic, start = y_end + gap / y_freq, frequency = y_freq)
  }
  
  n_tot <- k * p + m * (s + 1) + n
  
  pred_start <- stats::time(stats::ts(rep(NA, gap + 2), start = y_end, frequency = y_freq))[gap + 2]
  
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
      #
      # Start values, where given, take the place of the sample: the lags of
      # the first forecast period are their last p observations.
      temp_p <- stats::embed(as.matrix(if (is.null(y_start)) y else y_start), p)
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
      # A series that does not reach the forecast periods made stats::window()
      # warn and return fewer rows, which then failed to fit into 'x' with a
      # message that said nothing about what was missing.
      tsp_d <- stats::tsp(deterministic)
      tol <- 0.1 / y_freq
      if (abs(tsp_d[3] - y_freq) > 1e-8 || tsp_d[1] > x_time[1] + tol ||
          tsp_d[2] < x_time[length(x_time)] - tol) {
        stop("Argument 'deterministic' must cover the forecast periods from ",
             .ts_period_label(x_time[1], y_freq), " to ", .ts_period_label(x_time[length(x_time)], y_freq),
             " at the frequency of the model (", y_freq, "). It covers ",
             .ts_period_label(tsp_d[1], tsp_d[3]), " to ", .ts_period_label(tsp_d[2], tsp_d[3]),
             " at frequency ", tsp_d[3], ".", call. = FALSE)
      }
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
  if (!is.null(y_start)) {
    result[["start"]] <- as.numeric(pred_start)
    result[["y"]] <- y_start
  }
  
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

# The start values of a forecast, argument 'y_start' of prepare_forecast_input(),
# checked against the model. A list with 'y', the values the forecast reads and
# the ones after the end of the sample, one column per endogenous variable in
# the order of the model, and 'gap', the number of periods from the end of the
# sample to the last start value. The forecast starts in the period after that.
.check_forecast_start_values <- function(object, y_start) {

  if (!stats::is.ts(y_start)) {
    stop("Argument 'y_start' must be a time-series object (class 'ts') with the endogenous ",
         "variables of the model, such as the result of stats::ts() or stats::window().",
         call. = FALSE)
  }

  # The samplers that complete a panel not observed whole start each draw's
  # forecast from its own completion and overwrite the lags written here.
  if (!is.null(object[["data"]][["train"]][["constraints"]])) {
    stop("Argument 'y_start' is not available for a model estimated with missing = \"estimate\": ",
         "its forecast starts from each draw's completed panel, which takes the place of any ",
         "start values. Add the periods to the data the model is estimated with instead, with ",
         "NA where a variable was not observed, and the forecast starts after them.",
         call. = FALSE)
  }

  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  endogen <- object[["model"]][["endogen"]]
  tsp_y <- stats::tsp(object[["data"]][["train"]][["y"]])
  y_freq <- tsp_y[3]
  tsp_s <- stats::tsp(y_start)

  if (abs(tsp_s[3] - y_freq) > 1e-8) {
    stop("Argument 'y_start' has ", frequency_name(tsp_s[3]), " data, but the model is estimated ",
         "on ", frequency_name(y_freq), " data.", call. = FALSE)
  }

  # By name, or by position where the series has no names and the right number
  # of columns
  values <- as.matrix(y_start)
  names_s <- colnames(values)
  if (!is.null(names_s) && all(endogen %in% names_s)) {
    values <- values[, endogen, drop = FALSE]
  } else if (is.null(names_s) && ncol(values) == k) {
    colnames(values) <- endogen
  } else {
    stop("Argument 'y_start' must contain the endogenous variables of the model, ",
         paste0("'", endogen, "'", collapse = ", "), ", as columns of that name. It has ",
         if (is.null(names_s)) paste(ncol(values), "unnamed column(s)") else
           paste0("columns ", paste0("'", names_s, "'", collapse = ", ")), ".", call. = FALSE)
  }

  gap <- (tsp_s[2] - tsp_y[2]) * y_freq
  if (abs(gap - round(gap)) > 0.1) {
    stop("The periods of argument 'y_start' do not line up with those of the estimation sample.",
         call. = FALSE)
  }
  gap <- as.integer(round(gap))
  if (gap < 0) {
    stop("Argument 'y_start' ends in ", .ts_period_label(tsp_s[2], y_freq), ", before the end of ",
         "the estimation sample in ", .ts_period_label(tsp_y[2], y_freq), ". A forecast starts in ",
         "the last period of the sample or after it. A forecast from inside the sample is that of ",
         "a model estimated on a shorter one: create it from stats::window(data, end = ...).",
         call. = FALSE)
  }

  # The rows the forecast reads: the last p, the lags of its first period. The
  # last one at least, which is where the forecast starts in time, when there
  # are no lags.
  n_used <- max(p, 1)
  tt_s <- nrow(values)
  if (tt_s < n_used) {
    stop("Argument 'y_start' must contain the last ", n_used, " period(s) before the first ",
         "forecast period, ", .ts_period_label(tsp_s[2] - (n_used - 1) / y_freq, y_freq), " to ",
         .ts_period_label(tsp_s[2], y_freq), ", whose values are its lagged regressors (p = ", p,
         "). It has ", tt_s, ".", call. = FALSE)
  }
  if (p > 0) {
    used <- values[tt_s - p + 1:p, , drop = FALSE]
    if (anyNA(used)) {
      stop("Argument 'y_start' has missing values in its last ", p, " period(s), ",
           .ts_period_label(tsp_s[2] - (p - 1) / y_freq, y_freq), " to ",
           .ts_period_label(tsp_s[2], y_freq), ", in ",
           paste0("'", colnames(used)[colSums(is.na(used)) > 0], "'", collapse = ", "),
           ". They are the lags the forecast starts from and cannot be missing: fill them, ",
           "with a nowcast for example, or end 'y_start' where all variables are observed.",
           call. = FALSE)
    }
  }

  # Kept as well: the periods after the end of the sample, which predict() shows
  # as the history the forecast continues
  n_keep <- min(tt_s, max(n_used, gap))
  y <- stats::ts(values[tt_s - n_keep + 1:n_keep, , drop = FALSE], end = tsp_s[2],
                 frequency = y_freq)

  list("y" = y, "gap" = gap)
}

# Seasonal dummies of the estimation sample, 'season' with one row per period,
# continued over the last period of the sample and the 'n_ahead' after it. The
# pattern repeats itself with a period, the shortest shift under which the
# sample equals itself, and a period after the sample takes the values of the
# last period of the sample in the same season. A pattern without a period
# falls back on the first row that equals the last one and the rows after it.
.continue_seasonal_pattern <- function(season, n_ahead) {

  season <- as.matrix(season)
  tt <- nrow(season)
  j <- 0:n_ahead

  period <- NULL
  if (tt > 1) {
    for (l in 1:(tt - 1)) {
      if (all(season[(l + 1):tt, , drop = FALSE] == season[1:(tt - l), , drop = FALSE])) {
        period <- l
        break
      }
    }
  }
  if (!is.null(period)) {
    return(season[tt + j - period * ceiling(j / period), , drop = FALSE])
  }

  result <- matrix(NA_real_, length(j), ncol(season))
  last_obs <- season[tt, ]
  for (i in 1:tt) {
    if (identical(last_obs, season[i, ])) {
      rows <- i + j
      result[rows <= tt, ] <- season[rows[rows <= tt], , drop = FALSE]
      break
    }
  }
  result
}

# The number of periods between the end of the estimation sample and the period
# before the first forecast period: zero for a forecast that starts right after
# the sample, which is all of them unless add_forecast_input() was given
# 'y_start'.
.forecast_gap <- function(object) {
  start <- object[["model"]][["forecast_start"]]
  tsp_y <- stats::tsp(object[["data"]][["train"]][["y"]])
  if (is.null(start) || is.null(tsp_y)) {
    return(0L)
  }
  as.integer(round((start - tsp_y[2]) * tsp_y[3])) - 1L
}

# Refuses 'what' for a forecast that does not start right after the estimation
# sample, saying why.
.refuse_shifted_forecast <- function(object, what, why) {
  if (.forecast_gap(object) > 0) {
    stop("The forecast of this model starts after the end of its estimation sample, from the ",
         "values given in argument 'y_start' of add_forecast_input(), and ", what, " is not ",
         "available for it: ", why, call. = FALSE)
  }
  invisible(NULL)
}
