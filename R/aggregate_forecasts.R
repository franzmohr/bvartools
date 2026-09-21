#' Aggregate Forecasts to Annual Figures
#'
#' Turns the forecast draws of models estimated on quarterly or monthly data into
#' draws of the annual figures they imply, so that the models can be compared with
#' annual forecasts, such as the projections of international institutions, on the
#' same annual basis.
#'
#' @param object an object of class 'bvarmodel', 'bvecmodel', 'expandingwindow' or
#' 'modellist', whose models contain forecasts, i.e. \code{\link{add_posterior_forecasts}}
#' was already applied to them.
#' @param type a named character vector, which specifies how each variable is
#' aggregated. The names are endogenous variables of the models and the elements
#' one of \code{"growth"}, \code{"loglevel"} or \code{"level"}. See 'Details'.
#' Endogenous variables, which are not named, are dropped from the comparison.
#' @param scale the factor, by which the logarithms of variables of type
#' \code{"growth"} and \code{"loglevel"} are multiplied in the data. The default
#' of 100 corresponds to log levels and log differences in percent.
#'
#' @details
#' The three types describe how a variable is measured in the models and determine
#' how its annual figure is obtained:
#' \describe{
#'   \item{\code{"growth"}}{the change of a log level from one period to the next,
#'   such as \code{100 * diff(log(gdp))}. The changes are cumulated to log levels.}
#'   \item{\code{"loglevel"}}{a log level, such as \code{100 * log(gdp)} in a VEC
#'   model.}
#'   \item{\code{"level"}}{a variable, whose annual figure is its average over the
#'   year, such as an unemployment rate or an interest rate.}
#' }
#' For the first two the annual figure is the growth rate in percent of the annual
#' average of the levels, \code{100 * (mean(exp(cur / scale)) / mean(exp(prev / scale)) - 1)},
#' where \code{cur} and \code{prev} are the log levels of the periods of the year
#' and of the year before. This is also the growth rate of the annual sum of the
#' levels, since both years have the same number of periods, so it serves flows
#' such as GDP, whose annual figure is a sum, and prices, whose annual figure is an
#' average, alike. For \code{"level"} the annual figure is the average of the
#' periods of the year.
#'
#' The periods of a year, which were already observed at the end of the training
#' sample of a model, are taken from its data, and the remaining periods from each
#' draw of its forecast, so that every draw of the forecast becomes a draw of the
#' annual figure. Horizon 1 is the year of the forecast origin, i.e. the year of the
#' first period after the training sample, and horizon 2 the year after it. The
#' number of annual horizons is the number of years that the forecast horizon
#' covers from every origin, \code{floor(n_ahead / frequency)}, so that a
#' quarterly model, which forecasts eight quarters, provides the current and the
#' next year.
#'
#' The data of the models are aggregated in the same way and replace their data,
#' so that the realised annual values, against which the aggregated forecasts are
#' scored, are in \code{data$test$y}. \code{\link{add_forecast_errors}} uses them,
#' if its argument \code{test_sample} is omitted. Otherwise, \code{test_sample}
#' must be an annual time series, for example of official annual figures.
#'
#' The result can be passed as argument \code{object} to
#' \code{\link{create_external_forecast}}, which then reads the periods of the
#' external forecasts as years, but matches their publications to the training
#' samples of the models at the frequency of the data. The external forecasts are
#' scored against the same realised annual values as the models.
#'
#' The aggregated models only contain the forecasts and the data, which are
#' required to evaluate them. Therefore, aggregation is the last step before
#' \code{\link{add_forecast_errors}} and \code{\link{selection_criteria}}, which
#' only provides out-of-sample statistics for them. Aggregated VEC models become
#' objects of class 'bvarmodel', because their forecasts are those of the levels.
#'
#' @return An object of the same class as \code{object}, whose models contain the
#' draws of the annual forecasts in \code{posterior$forecast$forecasts}.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' # Create model
#' model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "const",
#'                           iterations = 10, burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#'
#' model <- use_expanding_window(model, start = 2005)
#' model <- add_priors(model, coef = list(v_i = 0.1, v_i_det = 0.01),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#' model <- add_forecast_input(model, n_ahead = 8)
#' model <- add_posterior_forecasts(model)
#'
#' # Annual averages of inflation and the interest rate
#' annual <- aggregate_forecasts(model, type = c(Dp = "level", r = "level"))
#'
#' # Artificial annual forecasts of the current and the next year, published in
#' # the middle of each quarter
#' fcst <- expand.grid(origin = 2005 + (0:7) / 4 + 0.1, year = 0:1,
#'                     variable = c("Dp", "r"), stringsAsFactors = FALSE)
#' fcst[["period"]] <- floor(fcst[["origin"]]) + fcst[["year"]]
#' fcst[["value"]] <- 2
#'
#' ext <- create_external_forecast(fcst, annual, data_lag = 1)
#'
#' # Both are scored against the annual averages of the data
#' models <- add_forecast_errors(combine_models(annual, ext))
#' selection_criteria(models)
#'
#' @family model comparison
#' @export
aggregate_forecasts <- function(object, type, scale = 100) {

  if (is.list(type)) {
    type <- unlist(type)
  }
  if (!is.character(type) || length(type) == 0 || is.null(names(type)) ||
      any(is.na(names(type)) | names(type) == "") || anyDuplicated(names(type)) > 0) {
    stop("Argument 'type' must be a named character vector with one element per variable,\nsuch as c(dy = \"growth\", u = \"level\").")
  }
  unknown <- !type %in% c("growth", "loglevel", "level")
  if (any(unknown)) {
    stop("The elements of argument 'type' must be 'growth', 'loglevel' or 'level', not ",
         paste0("'", unique(type[unknown]), "'", collapse = ", "), ".")
  }
  if (length(scale) != 1 || !is.numeric(scale) || !is.finite(scale) || scale <= 0) {
    stop("Argument 'scale' must be a single positive number.")
  }

  return(.aggregate_forecasts(object, type, scale))
}

# Walks a model object down to its models, the windows of an expanding window
# and the elements of a list of models, and aggregates each.
.aggregate_forecasts <- function(object, type, scale) {

  if ("externalforecast" %in% class(object)) {
    stop("Argument 'object' contains external forecasts, which are not aggregated. Aggregate\nthe forecasts of the models and pass the annual external forecasts to\n'create_external_forecast' together with them.")
  }

  if (any(c("modellist", "expandingwindow") %in% class(object))) {
    for (i in seq_along(object)) {
      object[[i]] <- .aggregate_forecasts(object[[i]], type, scale)
    }
    return(object)
  }

  if (!any(c("bvarmodel", "bvecmodel") %in% class(object))) {
    stop("Argument 'object' must be an object of class 'bvarmodel', 'bvecmodel',\n'expandingwindow' or 'modellist'.")
  }

  return(.aggregate_window(object, type, scale))
}

# The annual forecasts of one model.
.aggregate_window <- function(object, type, scale) {

  if (!is.null(object[["model"]][["aggregation"]])) {
    stop("The forecasts of the models in argument 'object' are already aggregated.")
  }

  draws <- .forecast_draws(object)
  if (is.null(draws)) {
    stop("At least one model in argument 'object' does not contain forecasts. Use\n'add_posterior_forecasts' before this function.")
  }

  # The forecasts of a VEC model are those of the levels, so the data of its
  # VAR representation are the ones they continue
  source <- object
  if ("bvecmodel" %in% class(object)) {
    source <- .vec_level_form(object)
  }

  y <- source[["data"]][["train"]][["y"]]
  tsp_y <- stats::tsp(y)
  freq <- tsp_y[3]
  if (freq <= 1 || abs(freq - round(freq)) > 1e-8) {
    stop("Forecasts can only be aggregated to annual figures, if the data have a whole\nnumber of periods per year above one, such as quarterly or monthly data.")
  }
  freq <- round(freq)

  endogen <- object[["model"]][["endogen"]]
  if (is.null(endogen)) {
    endogen <- dimnames(y)[[2]]
  }
  missing_vars <- !names(type) %in% endogen
  if (any(missing_vars)) {
    stop("Argument 'type' names the variable(s) ",
         paste0("'", names(type)[missing_vars], "'", collapse = ", "),
         ", which are not endogenous variables of the models in argument 'object'.")
  }

  k <- length(endogen)
  h <- object[["model"]][["h"]]
  if (is.null(h)) {
    h <- ncol(draws) %/% k
  }
  n_years <- h %/% freq
  if (n_years < 1) {
    stop("A forecast horizon of ", h, " periods does not cover a whole year of ",
         frequency_name(freq), " data from every forecast origin. Use\n'add_forecast_input' with 'n_ahead' of at least ", freq, " periods.")
  }

  # Periods are counted as whole periods of the data, so that they can be
  # compared exactly
  end <- round(tsp_y[2] * freq)
  first_year <- (end + 1) %/% freq
  years <- first_year + 0:(n_years - 1)

  original <- source[["data"]][["original"]][["endogen"]]
  if (is.null(original)) {
    original <- y
  }
  original <- .as_named_mts(original, endogen)
  period <- round(stats::time(original) * freq)

  # The year before the forecast origin is required for growth rates, from
  # its first period on
  observed <- period <= end & period >= (first_year - 1) * freq 
  n_draws <- nrow(draws)
  result <- matrix(NA_real_, n_draws, n_years * length(type))
  for (j in seq_along(type)) {
    var <- names(type)[j]
    pos <- which(endogen == var)
    past <- matrix(as.numeric(original[observed, var]), sum(observed), n_draws)
    future <- t(matrix(as.numeric(draws[, (0:(h - 1)) * k + pos]), n_draws, h))
    annual <- .annual_values(rbind(past, future), c(period[observed], end + 1:h),
                             years, type[j], scale, freq)
    if (anyNA(annual)) {
      stop("The data of variable '", var, "' do not cover the periods of the year ",
           "before the\nforecast origin ", .format_period(end + 1, freq), ", which the annual figures require.")
    }
    result[, (0:(n_years - 1)) * length(type) + j] <- t(annual)
  }
  dimnames(result) <- list(NULL, paste0(rep(names(type), n_years), "_",
                                        rep(1:n_years, each = length(type))))
  mc_stats <- coda::mcpar(draws)
  result <- coda::mcmc(result, start = mc_stats[1], end = mc_stats[2], thin = mc_stats[3])

  # The data, from which the realised values are taken, in the same annual terms
  annual_data <- .annual_data(original, type, scale, freq)

  object <- source
  object[["model"]][["k"]] <- length(type)
  object[["model"]][["endogen"]] <- names(type)
  object[["model"]][["h"]] <- n_years
  object[["model"]][["aggregation"]] <- list(type = type, scale = scale,
                                             frequency = freq, end = tsp_y[2])
  object[["data"]] <- list(original = list(endogen = annual_data),
                           train = list(y = stats::window(annual_data, end = first_year - 1)),
                           test = list(y = .realised_years(annual_data, years)))
  object[["posterior"]] <- list(forecast = list(forecasts = result))
  if (is.null(object[["data"]][["test"]][["y"]])) {
    object[["data"]][["test"]] <- NULL
  }
  class(object) <- c("bvarmodel", "list")

  return(object)
}

# Annual figures from the values of consecutive periods, one column per draw.
#
# values: periods x draws, the variable as the models see it
# period: the whole period of each row
# years:  the years whose figures are required
#
# Returns years x draws, NA where a period that a figure requires is missing.
.annual_values <- function(values, period, years, type, scale, freq) {

  result <- matrix(NA_real_, length(years), ncol(values))
  for (i in seq_along(years)) {
    current <- years[i] * freq + 0:(freq - 1)
    if (type == "level") {
      rows <- match(current, period)
      if (anyNA(rows)) {
        next
      }
      result[i, ] <- colMeans(values[rows, , drop = FALSE])
      next
    }

    if (type == "growth") {
      # The changes from the first period of the previous year on, cumulated
      # onto an arbitrary log level of that period, which cancels in the ratio
      rows <- match(current[1] - freq + 1:(2 * freq - 1), period)
      if (anyNA(rows)) {
        next
      }
      levels <- rbind(0, apply(values[rows, , drop = FALSE], 2, cumsum))
    } else {
      rows <- match(current[1] - freq + 0:(2 * freq - 1), period)
      if (anyNA(rows)) {
        next
      }
      # Relative to the first period, which keeps exp() within range
      levels <- sweep(values[rows, , drop = FALSE], 2, values[rows[1], ])
    }
    levels <- exp(levels / scale)
    result[i, ] <- 100 * (colMeans(levels[freq + 1:freq, , drop = FALSE]) /
                            colMeans(levels[1:freq, , drop = FALSE]) - 1)
  }

  return(result)
}

# The annual figures of the data, one column per variable, from the year before
# the first period of the data to the last year the data cover completely.
.annual_data <- function(original, type, scale, freq) {

  period <- round(stats::time(original) * freq)
  years <- (period[1] %/% freq - 1):((period[length(period)] + 1) %/% freq - 1)

  result <- matrix(NA_real_, length(years), length(type))
  for (j in seq_along(type)) {
    values <- as.numeric(original[, names(type)[j]])
    avail <- !is.na(values)
    result[, j] <- .annual_values(matrix(values[avail], ncol = 1), period[avail],
                                  years, type[j], scale, freq)
  }
  dimnames(result) <- list(NULL, names(type))

  return(stats::ts(result, start = years[1], frequency = 1))
}

# The realised values of the years of a forecast, one row per year. Only the
# years the data cover, which are the first ones, and NULL if there are none.
.realised_years <- function(annual_data, years) {
  pos <- match(years, round(stats::time(annual_data)))
  pos <- pos[!is.na(pos)]
  if (length(pos) == 0) {
    return(NULL)
  }
  return(as.matrix(annual_data)[pos, , drop = FALSE])
}

# A time series as a matrix with column names, also where it has one column.
.as_named_mts <- function(x, names) {
  tsp_x <- stats::tsp(x)
  x <- as.matrix(x)
  if (is.null(dimnames(x)[[2]])) {
    dimnames(x) <- list(NULL, names)
  }
  return(stats::ts(x, start = tsp_x[1], frequency = tsp_x[3]))
}

# A whole period of the data as, e.g., 2020Q3.
.format_period <- function(period, freq) {
  suffix <- switch(as.character(freq), "4" = "Q", "12" = "M", "2" = "H", ".")
  paste0(period %/% freq, suffix, period %% freq + 1)
}
