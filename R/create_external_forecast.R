#' Objects for Externally Produced Forecasts
#'
#' Turns externally produced point forecasts into an object, which can be evaluated
#' and compared with the forecasts of Bayesian VAR and VEC models.
#'
#' @param forecasts a data frame in long format, which contains the external point
#' forecasts. See 'Details' for the required columns.
#' @param object for \code{create_external_forecast} a model object of class
#' 'bvarmodel', 'bvecmodel', 'expandingwindow' or 'modellist', which is used as the
#' reference of the comparison. It provides the endogenous variables, the frequency
#' of the data and the ends of the training samples, to which the external forecasts
#' are matched. For the methods of the estimation workflow an object of class
#' 'externalforecast'.
#' @param n_ahead the maximum forecast horizon that is considered. If \code{NULL}
#' (default), the forecast horizon of the models in \code{object} is used, which
#' requires that \code{\link{add_forecast_input}} was already applied to them. For
#' models whose forecasts were aggregated with \code{\link{aggregate_forecasts}}
#' it is a number of years.
#' @param period name of the column of \code{forecasts}, which contains the period,
#' for which a forecast was made.
#' @param origin name of the column of \code{forecasts}, which contains the period,
#' in which a forecast was published.
#' @param variable name of the column of \code{forecasts}, which contains the names
#' of the forecasted variables. They must correspond to the names of the endogenous
#' variables of the models in \code{object}.
#' @param value name of the column of \code{forecasts}, which contains the values of
#' the forecasts.
#' @param by name of an optional column of \code{forecasts}, which contains the names
#' of the forecasters. If specified, one object is produced per forecaster and the
#' result is a list of class 'modellist'.
#' @param data_lag an integer specifying the number of periods, by which the
#' publication of the data of the endogenous variables lags behind. See 'Details'.
#' @param select either \code{"last"} (default) or \code{"first"} specifying which
#' forecast is used, if multiple publications are matched to the same training
#' sample. See 'Details'.
#'
#' @details
#' Argument \code{forecasts} must be a data frame in long format, where each row
#' contains a single point forecast. The names of the required columns can be
#' specified in the arguments \code{period}, \code{origin}, \code{variable} and
#' \code{value}. Periods can be provided either as objects of class 'Date' or as
#' numerics, which follow the convention of \code{\link[stats]{time}}, i.e. 2007.25
#' for the second quarter of 2007. The periods, for which a forecast was made, are
#' rounded to the frequency of the data of the models in \code{object}.
#'
#' Rounding places a forecast in the period of the data it falls into; it does
#' not convert between frequencies. Forecasts made at a frequency other than that
#' of the data are therefore refused: annual forecasts of a quarterly model,
#' whose periods such as 2020 and 2021 would otherwise be read as the first
#' quarters of those years and an annual growth rate scored as a quarterly one,
#' or monthly forecasts of a quarterly model, of which three would fall into the
#' same quarter. The frequency of the forecasts is read off the spacing of the
#' periods within a publication and, where every publication forecasts a single
#' period, off the position of the periods within the year.
#'
#' Annual forecasts are compared with models estimated on quarterly or monthly
#' data by aggregating the forecasts of the models to annual figures with
#' \code{\link{aggregate_forecasts}} first and passing the result as argument
#' \code{object}. The periods of the forecasts are then read as years, while
#' publications are still matched to the training samples at the frequency of the
#' data, and \code{data_lag} still counts periods of the data. Horizon 1 is the
#' year of the first period after the training sample, so that a forecast for the
#' year of its publication is a forecast of horizon 1 and one for the year after
#' it of horizon 2, whichever quarter it was published in. The realised values,
#' against which the forecasts are scored, are the annual figures of the data of
#' the models, which the aggregated forecasts of the models are scored against,
#' too, so that \code{\link{add_forecast_errors}} does not need argument
#' \code{test_sample}.
#'
#' In contrast to a model, an external forecast does not have a training sample.
#' Therefore, each publication is matched to the training sample, which ends closest
#' before the publication of the forecast, so that a model and an external forecaster
#' are evaluated on comparable information sets. Since the data of the endogenous
#' variables are usually published with a delay, argument \code{data_lag} can be used
#' to specify the number of periods, by which the last available observation lags
#' behind the publication of a forecast. Thus, a publication is matched to the last
#' training sample, which does not end after \code{origin - data_lag} periods. With
#' the default of one period an annual forecast, which was published in the course of
#' 2021, is matched to a training sample that ends in 2020, so that the forecast for
#' 2021 is a one-step ahead forecast.
#'
#' If multiple publications are matched to the same training sample, only one of them
#' is used, because otherwise a forecaster would enter the comparison multiple times
#' for the same period. Argument \code{select} controls whether the latest
#' (\code{"last"}) or the earliest (\code{"first"}) of those publications is used.
#'
#' The resulting object mimics the structure of a model, which was estimated with the
#' expanding window approach of \code{\link{use_expanding_window}}, where each
#' publication corresponds to one estimation window. It can be added to a list of
#' models with \code{\link{combine_models}} and the functions, which are applied to
#' obtain posterior draws, are without effect for it. Since external forecasts are
#' point forecasts, their forecast errors consist of a single draw. Accordingly, the
#' credible bands of the out-of-sample statistics in \code{\link{selection_criteria}}
#' are degenerate and in-sample criteria are not available.
#'
#' @return A list of class 'externalforecast' or, if argument \code{by} is specified,
#' a list of class 'modellist', which contains one object of class 'externalforecast'
#' per forecaster.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' # Create model
#' model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "none",
#'                           error = "gamma", iterations = 10, burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#'
#' model <- use_expanding_window(model, start = 2007)
#'
#' # Artificial external forecasts of two forecasters
#' fcst <- expand.grid(origin = c(2007, 2007.25),
#'                     h = 1:2,
#'                     variable = c("Dp", "r"),
#'                     forecaster = c("A", "B"),
#'                     stringsAsFactors = FALSE)
#' fcst[["period"]] <- fcst[["origin"]] + fcst[["h"]] / 4
#' fcst[["value"]] <- 0
#'
#' # Create objects of the external forecasts, where the data of the endogenous
#' # variables are assumed to be published with a delay of one quarter
#' ext <- create_external_forecast(fcst, model, n_ahead = 4,
#'                                 by = "forecaster", data_lag = 1)
#'
#' # Calculate forecast errors
#' ext <- add_forecast_errors(ext, test_sample = us_macrodata)
#'
#' # Compare the forecast performance
#' selection_criteria(ext)
#'
#' @family model comparison
#' @export
create_external_forecast <- function(forecasts, object, n_ahead = NULL,
                                     period = "period", origin = "origin",
                                     variable = "variable", value = "value",
                                     by = NULL, data_lag = 1, select = "last") {

  if (!is.data.frame(forecasts)) {
    stop("Argument 'forecasts' must be a data frame.")
  }

  required <- c(period, origin, variable, value, by)
  missing_cols <- required[!required %in% names(forecasts)]
  if (length(missing_cols) > 0) {
    stop("Argument 'forecasts' does not contain the column(s) ",
         paste0("'", missing_cols, "'", collapse = ", "), ".")
  }

  if (!select %in% c("last", "first")) {
    stop("Argument 'select' must be either 'last' or 'first'.")
  }

  if (length(data_lag) != 1 || !is.finite(data_lag)) {
    stop("Argument 'data_lag' must be a single finite number.")
  }

  # Specifications of the models, against which the external forecasts are compared
  ref <- get_reference_periods(object)

  if (is.null(n_ahead)) {
    n_ahead <- ref[["h"]]
    if (is.null(n_ahead)) {
      stop("Argument 'n_ahead' must be specified, because argument 'object' does not\ncontain a forecast horizon. You might want to use function 'add_forecast_input'\nbefore this function.")
    }
  }
  n_ahead <- as.integer(n_ahead)

  freq <- ref[["frequency"]]
  endogen <- ref[["endogen"]]
  # Annual, if the forecasts of the models were aggregated, and the frequency of
  # the data otherwise
  period_freq <- ref[["period_frequency"]]

  # The periods of the forecasts are rounded to the frequency of the forecasts,
  # whereas the publications are not, so that publications within a period can
  # be ordered
  fcst <- data.frame(period = to_model_time(forecasts[[period]], period_freq, round = TRUE),
                     origin = to_model_time(forecasts[[origin]], freq, round = FALSE),
                     variable = as.character(forecasts[[variable]]),
                     value = as.numeric(forecasts[[value]]),
                     stringsAsFactors = FALSE)

  if (is.null(by)) {
    groups <- rep("", nrow(fcst))
  } else {
    groups <- as.character(forecasts[[by]])
  }

  # Rows, which cannot be used for the comparison. Incomplete rows are omitted
  # silently, whereas variables, which are not part of the models, might indicate
  # that the wrong forecasts were provided
  unknown <- !is.na(fcst[, "variable"]) & !fcst[, "variable"] %in% endogen
  if (any(unknown)) {
    warning("Argument 'forecasts' contains the variable(s) ",
            paste0("'", unique(fcst[unknown, "variable"]), "'", collapse = ", "),
            ", which are not endogenous variables of the models in argument 'object'.\nThey are omitted.")
  }
  keep <- fcst[, "variable"] %in% endogen & !is.na(fcst[, "period"]) &
    !is.na(fcst[, "origin"]) & !is.na(fcst[, "value"])
  fcst <- fcst[keep, ]
  groups <- groups[keep]

  if (nrow(fcst) == 0) {
    stop("Argument 'forecasts' does not contain any usable forecast.")
  }

  check_forecast_frequency(fcst, groups, period_freq)

  group_names <- unique(groups)

  # Forecasters without a usable publication produce an empty element, which must
  # be maintained in the list so that it can be reported below
  result <- vector("list", length(group_names))
  names(result) <- group_names
  for (i in seq_along(group_names)) {
    result[i] <- list(build_external_forecast(fcst[groups == group_names[i], ],
                                              ref = ref, n_ahead = n_ahead,
                                              data_lag = data_lag, select = select))
  }

  empty <- unlist(lapply(result, is.null))
  if (all(empty)) {
    stop("None of the forecasts in argument 'forecasts' could be matched to a training\nsample of the models in argument 'object'. You might want to check arguments\n'origin' and 'data_lag'.")
  }
  if (any(empty)) {
    warning("No forecast of ", paste0("'", group_names[empty], "'", collapse = ", "),
            " could be matched to a training sample of the models in argument 'object'.")
    result <- result[!empty]
  }

  if (is.null(by)) {
    return(result[[1]])
  }

  class(result) <- c("modellist", "list")

  return(result)
}

#' Refuses Forecasts at Another Frequency Than the Data
#'
#' Rounding puts a forecast into the period of the data it falls into, which is
#' what lets a calendar date stand for a quarter. It is not a conversion between
#' frequencies, and a forecast at another frequency went through it without a
#' word: an annual forecast of 2021 became a forecast of 2021Q1 of a quarterly
#' model, and an annual growth rate was scored as a quarterly one.
#'
#' A single period cannot give that away -- 2021 is also a first quarter -- so
#' the frequency is read off the forecasts together. Within one publication of
#' one variable, periods that are further apart than one period of the data are
#' of a lower frequency, and two that round to the same period are of a higher
#' one. Where every publication forecasts a single period, periods that all sit
#' at the same position within the year, over more than one year, are of a lower
#' frequency too.
#'
#' @param fcst the data frame of \code{create_external_forecast}, with the
#' periods already rounded to the frequency of the data.
#' @param groups the forecaster of each row of \code{fcst}.
#' @param freq the frequency of the data of the models.
#'
#' @return \code{NULL}, invisibly. Called for the error.
#'
#' @noRd
check_forecast_frequency <- function(fcst, groups, freq) {

  # Whole periods of the data, so that the comparisons below are exact
  index <- round(fcst[, "period"] * freq)
  publication <- list(groups, fcst[, "origin"], fcst[, "variable"])

  steps <- unlist(lapply(split(index, publication, drop = TRUE), function(x) {
    diff(sort(x))
  }), use.names = FALSE)

  if (any(steps == 0)) {
    stop("Argument 'forecasts' contains several forecasts of the same variable in one ",
         "publication that fall into the same period of the data, which is ",
         frequency_name(freq), ". Forecasts at a higher frequency than the data ",
         "are not aggregated; supply them at the frequency of the data.")
  }

  if (length(steps) > 0) {
    step <- min(steps)
  } else if (freq > 1) {
    # One period per publication. Periods that never leave one position within
    # the year, over several years, are of a frequency of once a year.
    position <- index %% freq
    years <- unique(index %/% freq)
    step <- if (length(unique(position)) == 1 && length(years) > 1) freq else 1
  } else {
    step <- 1
  }

  if (step > 1) {
    # Annual forecasts can be compared with the annual figures that the
    # forecasts of the models imply
    remedy <- "."
    if (isTRUE(all.equal(freq / step, 1))) {
      remedy <- paste0(", or aggregate the forecasts of the models to annual figures ",
                       "with 'aggregate_forecasts' and pass the result as argument 'object'.")
    }
    stop("The periods in argument 'forecasts' are ", step, " periods of the data ",
         "apart and the data are ", frequency_name(freq), ", so they appear to be ",
         frequency_name(freq / step), " forecasts. Periods are only rounded to the ",
         "frequency of the data, not converted to it: a forecast of a year would be ",
         "compared with the first period of that year. Supply forecasts at the ",
         "frequency of the data", remedy)
  }

  invisible(NULL)
}

#' The Name of a Frequency
#'
#' @param freq a number of periods per year.
#'
#' @return A character string.
#'
#' @noRd
frequency_name <- function(freq) {
  names <- c("1" = "annual", "2" = "semi-annual", "4" = "quarterly", "12" = "monthly")
  name <- names[as.character(freq)]
  if (is.na(name)) {
    return(paste(format(freq), "periods per year"))
  }
  unname(name)
}

#' Reference Periods of a Model Object
#'
#' Obtains the ends of the training samples and further specifications, which are
#' required to match external forecasts to the models of a comparison.
#'
#' @param object an object of class 'bvarmodel', 'bvecmodel', 'expandingwindow' or
#' 'modellist'.
#'
#' @return A list with the ends of the training samples, the frequency of the data,
#' the frequency of the forecasts -- annual for models whose forecasts were
#' aggregated with \code{aggregate_forecasts}, the frequency of the data otherwise
#' --, the forecast horizon in periods of the forecasts, the names of the
#' endogenous variables and their data at the frequency of the forecasts.
#'
#' @noRd
get_reference_periods <- function(object) {

  if (any(c("modellist", "expandingwindow") %in% class(object))) {

    temp <- lapply(object, get_reference_periods)

    freq <- unique(unlist(lapply(temp, function(x) {x[["frequency"]]})))
    if (length(freq) > 1) {
      stop("The models in argument 'object' do not have the same frequency.")
    }
    period_freq <- unique(unlist(lapply(temp, function(x) {x[["period_frequency"]]})))
    if (length(period_freq) > 1) {
      stop("The forecasts of some models in argument 'object' are aggregated to annual\nfigures and those of others are not. Use 'aggregate_forecasts' on all of them.")
    }

    h <- unlist(lapply(temp, function(x) {x[["h"]]}))
    if (length(h) > 0) {
      h <- max(h)
    } else {
      h <- NULL
    }

    # The model with the longest training sample provides the data of the
    # endogenous variables
    ends <- lapply(temp, function(x) {x[["ends"]]})
    pos <- which.max(unlist(lapply(ends, max)))

    return(list(ends = sort(unique(unlist(ends))),
                frequency = freq,
                period_frequency = period_freq,
                h = h,
                endogen = temp[[pos]][["endogen"]],
                y = temp[[pos]][["y"]]))
  }

  if (!any(c("bvarmodel", "bvecmodel") %in% class(object))) {
    stop("Argument 'object' must be an object of class 'bvarmodel', 'bvecmodel',\n'expandingwindow' or 'modellist'.")
  }

  y <- object[["data"]][["train"]][["y"]]
  if (is.null(y)) {
    stop("At least one model in argument 'object' does not contain a training sample.")
  }
  tsp_y <- stats::tsp(y)

  endogen <- object[["model"]][["endogen"]]
  if (is.null(endogen)) {
    endogen <- dimnames(y)[[2]]
  }

  # The original data are used for the training samples of the external forecasts,
  # because they can start before the training sample of a model
  orig <- object[["data"]][["original"]][["endogen"]]
  if (is.null(orig)) {
    orig <- y
  }

  # A model whose forecasts were aggregated to annual figures is matched to
  # publications by the end of the data it was estimated on, while its data are
  # the annual figures
  aggregation <- object[["model"]][["aggregation"]]
  if (!is.null(aggregation)) {
    return(list(ends = aggregation[["end"]],
                frequency = aggregation[["frequency"]],
                period_frequency = tsp_y[3],
                h = object[["model"]][["h"]],
                endogen = endogen,
                y = orig))
  }

  return(list(ends = tsp_y[2],
              frequency = tsp_y[3],
              period_frequency = tsp_y[3],
              h = object[["model"]][["h"]],
              endogen = endogen,
              y = orig))
}

#' Time of a Period in the Convention of Time-Series Objects
#'
#' Converts the periods of external forecasts into numerics, which follow the
#' convention of \code{\link[stats]{time}}.
#'
#' @param x a vector of periods, either of class 'Date', 'POSIXt', 'character' or
#' 'numeric'.
#' @param frequency the frequency of the data of the models of the comparison.
#' @param round logical. If \code{TRUE}, the result is rounded down to the frequency
#' of the data.
#'
#' @return A numeric vector.
#'
#' @noRd
to_model_time <- function(x, frequency, round = TRUE) {

  if (inherits(x, "POSIXt")) {
    x <- as.Date(x)
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    # Periods, which are provided as calendar dates, must not be interpreted as
    # numerics, because their format differs from the one of time-series objects.
    # Missing periods are maintained and omitted at a later stage
    avail <- !is.na(x)
    if (all(grepl("^[0-9]{4}-[0-9]{2}", x[avail]))) {
      temp <- rep(NA_character_, length(x))
      temp[avail] <- paste0(substring(x[avail], 1, 7), "-01")
      x <- as.Date(temp)
    } else {
      x <- suppressWarnings(as.numeric(x))
      if (any(is.na(x) & avail)) {
        stop("At least one period in argument 'forecasts' could neither be interpreted as a\ndate nor as a numeric.")
      }
    }
  }

  if (inherits(x, "Date")) {

    result <- rep(NA_real_, length(x))
    avail <- !is.na(x)
    x <- x[avail]
    year <- as.numeric(format(x, "%Y"))

    if (round) {
      month <- as.numeric(format(x, "%m"))
      result[avail] <- year + floor((month - 1) * frequency / 12) / frequency
    } else if (frequency >= 1 && abs(12 / frequency - round(12 / frequency)) < 1e-8) {
      # A date is placed inside the period it falls in, by the share of that
      # period's days that have passed. The share of the year's days would put it
      # in the wrong period near the boundaries, since periods of months are not
      # equally long: 1 April is day 91 of 365, a share of 0.247 of the year, and
      # would count as the first quarter.
      month <- as.numeric(format(x, "%m"))
      months <- round(12 / frequency)
      period <- floor((month - 1) / months)
      first <- year * 12 + period * months
      start <- as.Date(sprintf("%04d-%02d-01", first %/% 12, first %% 12 + 1))
      after <- first + months
      end <- as.Date(sprintf("%04d-%02d-01", after %/% 12, after %% 12 + 1))
      share <- as.numeric(x - start) / as.numeric(end - start)
      result[avail] <- year + (period + share) / frequency
    } else {
      # Frequencies that do not divide the year into months: publications are
      # ordered by their day of the year
      day <- as.numeric(format(x, "%j")) - 1
      len <- as.numeric(format(as.Date(paste0(year, "-12-31")), "%j"))
      result[avail] <- year + day / len
    }

    return(result)
  }

  x <- as.numeric(x)
  if (round) {
    x <- floor(x) + floor((x - floor(x)) * frequency + 1e-8) / frequency
  }

  return(x)
}

#' Object of the External Forecasts of a Single Forecaster
#'
#' @param fcst a data frame with the columns \code{"period"}, \code{"origin"},
#' \code{"variable"} and \code{"value"}, where the periods are already converted
#' into the convention of time-series objects.
#' @param ref the result of a call to \code{get_reference_periods}.
#' @param n_ahead the maximum forecast horizon.
#' @param data_lag the publication lag of the data of the endogenous variables.
#' @param select either \code{"last"} or \code{"first"}.
#'
#' @return A list of class 'externalforecast' or \code{NULL}, if no publication
#' could be matched to a training sample.
#'
#' @noRd
build_external_forecast <- function(fcst, ref, n_ahead, data_lag, select) {

  ends <- ref[["ends"]]
  freq <- ref[["frequency"]]
  period_freq <- ref[["period_frequency"]]
  aggregated <- period_freq != freq
  endogen <- ref[["endogen"]]
  k <- length(endogen)
  tol <- 1e-8

  origins <- sort(unique(fcst[, "origin"]))
  # The last observation, which was available at the time of a publication
  cutoff <- origins - data_lag / freq
  # Position of the training sample, which ends closest before a publication
  pos <- findInterval(cutoff + tol, ends)
  origins <- origins[pos > 0]
  pos <- pos[pos > 0]

  if (length(pos) == 0) {
    return(NULL)
  }

  # Only one publication per training sample is used, because otherwise a
  # forecaster would enter the comparison multiple times for the same period
  windows <- unique(pos)
  sel <- rep(NA_real_, length(windows))
  for (i in seq_along(windows)) {
    temp <- origins[pos == windows[i]]
    sel[i] <- ifelse(select == "last", max(temp), min(temp))
  }

  result <- list()
  for (i in seq_along(windows)) {

    end_i <- ends[windows[i]]
    temp <- fcst[fcst[, "origin"] == sel[i], ]

    # The last period of the forecasts that the data cover completely. For
    # annual forecasts and quarterly data that is the last complete year, so
    # that horizon 1 is the year of the forecast origin.
    end_period <- end_i
    if (aggregated) {
      end_period <- (round(end_i * freq) + 1) %/% freq - 1
    }

    # Forecast horizon of the periods of the publication
    h_i <- (temp[, "period"] - end_period) * period_freq
    valid <- abs(h_i - round(h_i)) < 1e-6
    h_i <- round(h_i)
    valid <- valid & h_i >= 1 & h_i <= n_ahead
    temp <- temp[valid, ]
    h_i <- h_i[valid]

    if (nrow(temp) == 0) {
      next
    }

    # The columns of the forecasts are ordered by the endogenous variables within
    # the forecast horizons, where variables without a forecast remain missing
    fcst_i <- matrix(NA_real_, 1, n_ahead * k)
    dimnames(fcst_i) <- list(NULL, paste0(rep(endogen, n_ahead), "_",
                                          rep(1:n_ahead, each = k)))
    fcst_i[1, (h_i - 1) * k + match(temp[, "variable"], endogen)] <- temp[, "value"]

    model_i <- list(type = "External",
                    algorithm = "external",
                    k = k,
                    p = 0L,
                    m = 0L,
                    s = 0L,
                    n = 0L,
                    varsel = "none",
                    endogen = endogen,
                    structural = FALSE,
                    tvp = FALSE,
                    h = n_ahead,
                    origin = sel[i])

    result_i <- list("model" = model_i,
                     "data" = list("original" = list("endogen" = ref[["y"]]),
                                   "train" = list("y" = stats::window(ref[["y"]], end = end_period))),
                     "posterior" = list("forecast" = list("forecasts" = coda::mcmc(fcst_i))))

    # Annual forecasts are scored against the annual figures of the data, as
    # the aggregated forecasts of the models are
    if (aggregated) {
      result_i[["model"]][["aggregation"]] <- list(frequency = freq, end = end_i)
      realised <- .realised_years(ref[["y"]], end_period + 1:n_ahead)
      if (!is.null(realised)) {
        result_i[["data"]][["test"]] <- list("y" = realised)
      }
    }

    class(result_i) <- c("externalwindow", "bvarmodel", "list")

    result <- c(result, list(result_i))
  }

  if (length(result) == 0) {
    return(NULL)
  }

  # The object mimics a model, which was estimated with the expanding window
  # approach, so that the existing methods of that class can be used
  class(result) <- c("externalforecast", "expandingwindow", "list")

  return(result)
}
