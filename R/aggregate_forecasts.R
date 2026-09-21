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
#' @param code a named integer vector of the transformation codes of FRED-MD and
#' FRED-QD, \code{1:7}, which were applied to the variables of the models -- see
#' \code{\link{transform_variables}} and 'Details'. The names are endogenous
#' variables of the models. Endogenous variables, which are not named, are
#' dropped from the comparison.
#' @param target either \code{"average"} (default) or \code{"q4q4"}, the annual
#' figure that is compared. See 'Details'.
#' @param levels an optional time-series object with the untransformed series, to
#' which \code{code} was applied, i.e. the argument \code{x} of
#' \code{\link{transform_variables}}. Its columns are named after the variables.
#' Required for codes 2, 3, 6 and 7. See 'Details'.
#' @param scale the factor, by which the transformed series of codes 4 to 7 were
#' multiplied before they were passed to the models, such as 100 for log
#' differences in percent. The default of 1 corresponds to the result of
#' \code{\link{transform_variables}}.
#'
#' @details
#' Every draw of a forecast is turned back into a path of the untransformed
#' series by reversing the transformation of its code, and continued from the
#' periods, which were observed at the end of the training sample. The periods of
#' a year, which were already observed, are thus taken from the data and the
#' remaining periods from the draw, so that every draw of the forecast becomes a
#' draw of the annual figure.
#'
#' The codes determine the annual figure. The multiplicative codes 4 to 7 --
#' logarithms, their differences and differences of growth rates -- describe a
#' series such as GDP or a price index, whose annual figure is a growth rate in
#' percent. The additive codes 1 to 3 describe a series such as an unemployment
#' rate or an interest rate, whose annual figure is a level. Argument
#' \code{target} determines, which growth rate and which level:
#' \describe{
#'   \item{\code{"average"}}{the growth of the annual average of the levels over
#'   the annual average of the year before for codes 4 to 7, and the annual
#'   average for codes 1 to 3. This is the convention of, e.g., the World
#'   Economic Outlook of the IMF. The growth of the annual average is also the
#'   growth of the annual sum, since both years have the same number of periods,
#'   so it serves flows such as GDP and averages such as prices alike.}
#'   \item{\code{"q4q4"}}{the growth of the level of the last period of the year
#'   over the last period of the year before, i.e. the fourth quarter over the
#'   fourth quarter or December over December, for codes 4 to 7, and the level of
#'   the last period of the year for codes 1 to 3. This is the convention of,
#'   e.g., the Summary of Economic Projections of the Federal Reserve.}
#' }
#'
#' Reversing a difference requires the level, from which it starts. For codes 2,
#' 3, 6 and 7 it is taken from argument \code{levels}, which is therefore
#' required for them. Codes 1, 4 and 5 can be reversed from the data of the
#' models, because the unknown level of a logarithm cancels from a growth rate.
#' If \code{levels} is provided, it must reproduce the data of the models under
#' \code{code} and \code{scale}, which guards against a wrong \code{scale}, and
#' it also provides the realised annual figures.
#'
#' Horizon 1 is the year of the forecast origin, i.e. the year of the first
#' period after the training sample, and horizon 2 the year after it. The number
#' of annual horizons is the number of years that the forecast horizon covers
#' from every origin, \code{floor(n_ahead / frequency)}, so that a quarterly
#' model, which forecasts eight quarters, provides the current and the next year.
#'
#' The realised annual figures, against which the aggregated forecasts are
#' scored, are obtained from \code{levels} or, without it, from the data of the
#' models in the same way, and are put in \code{data$test$y}.
#' \code{\link{add_forecast_errors}} uses them, if its argument
#' \code{test_sample} is omitted. Otherwise, \code{test_sample} must be an annual
#' time series, for example of official annual figures.
#'
#' The result can be passed as argument \code{object} to
#' \code{\link{create_external_forecast}}, which then reads the periods of the
#' external forecasts as years, but matches their publications to the training
#' samples of the models at the frequency of the data. The external forecasts are
#' scored against the same realised annual figures as the models.
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
#' # Inflation and the interest rate enter the model as they are
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
#' annual <- aggregate_forecasts(model, code = c(Dp = 1, r = 1))
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
#' @references
#'
#' McCracken, M. W., & Ng, S. (2021). FRED-QD: A quarterly database for
#' macroeconomic research. \emph{Federal Reserve Bank of St. Louis Review,
#' 103}(1), 1--44.
#'
#' @family model comparison
#' @export
aggregate_forecasts <- function(object, code, target = "average", levels = NULL,
                                scale = 1) {

  if (is.list(code)) {
    code <- unlist(code)
  }
  if (!is.numeric(code) || length(code) == 0 || is.null(names(code)) ||
      any(is.na(names(code)) | names(code) == "") || anyDuplicated(names(code)) > 0) {
    stop("Argument 'code' must be a named numeric vector of transformation codes with one\nelement per variable, such as c(gdp = 5, unrate = 1).")
  }
  if (any(is.na(code)) || any(!code %in% 1:7)) {
    stop("Argument 'code' must only contain the transformation codes 1 to 7. See\n'transform_variables'.")
  }
  code <- stats::setNames(as.integer(code), names(code))

  if (length(target) != 1 || !target %in% c("average", "q4q4")) {
    stop("Argument 'target' must be either 'average' or 'q4q4'.")
  }

  if (length(scale) != 1 || !is.numeric(scale) || !is.finite(scale) || scale <= 0) {
    stop("Argument 'scale' must be a single positive number.")
  }

  if (is.null(levels)) {
    need <- code %in% c(2, 3, 6, 7)
    if (any(need)) {
      stop("Codes 2, 3, 6 and 7 are differences, which are reversed from the level they\nstart from. Argument 'levels' must provide the untransformed series of ",
           paste0("'", names(code)[need], "'", collapse = ", "), ".")
    }
  } else {
    if (!"ts" %in% class(levels)) {
      stop("Argument 'levels' must be an object of class 'ts'.")
    }
    levels <- .as_named_mts(levels, names(code)[1])
    missing_vars <- !names(code) %in% dimnames(levels)[[2]]
    if (any(missing_vars)) {
      stop("Argument 'levels' does not contain the variable(s) ",
           paste0("'", names(code)[missing_vars], "'", collapse = ", "), ".")
    }
  }

  spec <- list(code = code, target = target, levels = levels, scale = scale)

  return(.aggregate_forecasts(object, spec))
}

# Walks a model object down to its models, the windows of an expanding window
# and the elements of a list of models, and aggregates each.
.aggregate_forecasts <- function(object, spec) {

  if ("externalforecast" %in% class(object)) {
    stop("Argument 'object' contains external forecasts, which are not aggregated. Aggregate\nthe forecasts of the models and pass the annual external forecasts to\n'create_external_forecast' together with them.")
  }

  if (any(c("modellist", "expandingwindow") %in% class(object))) {
    for (i in seq_along(object)) {
      object[[i]] <- .aggregate_forecasts(object[[i]], spec)
    }
    return(object)
  }

  if (!any(c("bvarmodel", "bvecmodel") %in% class(object))) {
    stop("Argument 'object' must be an object of class 'bvarmodel', 'bvecmodel',\n'expandingwindow' or 'modellist'.")
  }

  return(.aggregate_window(object, spec))
}

# The annual forecasts of one model.
.aggregate_window <- function(object, spec) {

  code <- spec[["code"]]
  scale <- spec[["scale"]]
  levels <- spec[["levels"]]

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
  missing_vars <- !names(code) %in% endogen
  if (any(missing_vars)) {
    stop("Argument 'code' names the variable(s) ",
         paste0("'", names(code)[missing_vars], "'", collapse = ", "),
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
  # The year before the forecast origin, which the growth rates of its year
  # require, and the period before the end of the sample, from which a second
  # difference is reversed
  start <- min((first_year - 1) * freq, end - 1)

  original <- source[["data"]][["original"]][["endogen"]]
  if (is.null(original)) {
    original <- y
  }
  original <- .as_named_mts(original, endogen)

  if (!is.null(levels)) {
    .check_levels(levels, original, code, scale, freq)
  }

  n_draws <- nrow(draws)
  result <- matrix(NA_real_, n_draws, n_years * length(code))
  for (j in seq_along(code)) {
    var <- names(code)[j]
    pos <- which(endogen == var)
    past <- .past_levels(var, code[[j]], original, levels, start, end, scale, freq)
    if (is.null(past)) {
      stop("The data of variable '", var, "' do not cover the periods from ",
           .format_period(start, freq), " to ", .format_period(end, freq),
           ", which the\nannual figures of the forecast origin ",
           .format_period(end + 1, freq), " require.")
    }
    future <- t(matrix(as.numeric(draws[, (0:(h - 1)) * k + pos]), n_draws, h))
    if (code[[j]] >= 4) {
      future <- future / scale
    }
    path <- .continue_levels(code[[j]], past[["x"]], future, past[["log_offset"]])
    annual <- .annual_figures(path, start:(end + h), years, code[[j]] >= 4,
                              spec[["target"]], freq)
    result[, (0:(n_years - 1)) * length(code) + j] <- t(annual)
  }
  dimnames(result) <- list(NULL, paste0(rep(names(code), n_years), "_",
                                        rep(1:n_years, each = length(code))))
  mc_stats <- coda::mcpar(draws)
  result <- coda::mcmc(result, start = mc_stats[1], end = mc_stats[2], thin = mc_stats[3])

  # The realised values, in the same annual terms
  annual_data <- .annual_data(original, levels, code, spec[["target"]], scale, freq)

  object <- source
  object[["model"]][["k"]] <- length(code)
  object[["model"]][["endogen"]] <- names(code)
  object[["model"]][["h"]] <- n_years
  object[["model"]][["aggregation"]] <- list(code = code, target = spec[["target"]],
                                             scale = scale, frequency = freq,
                                             end = tsp_y[2])
  object[["data"]] <- list(original = list(endogen = annual_data),
                           train = list(y = stats::window(annual_data, end = first_year - 1)))
  realised <- .realised_years(annual_data, years)
  if (!is.null(realised)) {
    object[["data"]][["test"]] <- list(y = realised)
  }
  object[["posterior"]] <- list(forecast = list(forecasts = result))
  class(object) <- c("bvarmodel", "list")

  return(object)
}

# The untransformed series up to the end of a training sample, from which its
# forecasts are continued, over the periods start to end. Taken from the levels
# where they are given. Without them the codes 1, 4 and 5 are reversed from the
# data of the model: code 1 is the level itself, and codes 4 and 5 give the
# level up to a factor, which cancels from a growth rate, and which is chosen
# so that the last level is one. For code 4, 'log_offset' is the logarithm of
# that factor, by which the log levels of the forecasts are shifted in turn.
#
# Returns NULL, if a required period is missing.
.past_levels <- function(var, code, original, levels, start, end, scale, freq) {

  periods <- start:end
  if (!is.null(levels)) {
    x <- .values_at(levels[, var], periods, freq)
    if (anyNA(x)) {
      return(NULL)
    }
    return(list(x = x, log_offset = 0))
  }

  z <- .values_at(original[, var], periods, freq)
  if (code == 1) {
    x <- z
    log_offset <- 0
  } else if (code == 4) {
    log_offset <- z[length(z)] / scale
    x <- exp(z / scale - log_offset)
  } else {
    # Code 5: the first change is the one into the first period, which is not
    # needed and may be missing
    z[1] <- 0
    x <- exp(cumsum(z / scale) - sum(z / scale))
    log_offset <- 0
  }
  if (anyNA(x)) {
    return(NULL)
  }

  return(list(x = x, log_offset = log_offset))
}

# The values of a time series in the given whole periods, NA where it has none.
.values_at <- function(x, periods, freq) {
  as.numeric(x)[match(periods, round(stats::time(x) * freq))]
}

# The untransformed paths of the draws of a forecast: the past levels followed
# by the forecast of each draw, with the transformation of its code reversed.
#
# past:       the levels up to the end of the training sample
# z:          periods x draws, the transformed forecasts, divided by the scale
#             for codes 4 to 7
# log_offset: for code 4, the logarithm of the factor that 'past' was divided by
#
# Returns (length(past) + nrow(z)) x draws.
.continue_levels <- function(code, past, z, log_offset = 0) {

  n <- length(past)
  x_t <- past[n]
  x_s <- past[n - 1]

  cumulate <- function(m) {
    m <- matrix(m, nrow(z), ncol(z))
    for (i in seq_len(nrow(m))[-1]) {
      m[i, ] <- m[i, ] + m[i - 1, ]
    }
    m
  }

  path <- switch(as.character(code),
                 "1" = z,
                 "2" = x_t + cumulate(z),
                 "3" = x_t + cumulate(x_t - x_s + cumulate(z)),
                 "4" = exp(z - log_offset),
                 "5" = x_t * exp(cumulate(z)),
                 "6" = x_t * exp(cumulate(log(x_t / x_s) + cumulate(z))),
                 "7" = {
                   growth <- 1 + x_t / x_s - 1 + cumulate(z)
                   for (i in seq_len(nrow(growth))[-1]) {
                     growth[i, ] <- growth[i, ] * growth[i - 1, ]
                   }
                   x_t * growth
                 })

  return(rbind(matrix(past, n, ncol(z)), matrix(path, nrow(z), ncol(z))))
}

# Annual figures from the levels of consecutive periods, one column per draw.
#
# values:         periods x draws, the untransformed series
# period:         the whole period of each row
# years:          the years whose figures are required
# multiplicative: TRUE for a growth rate in percent, FALSE for a level
# target:         "average" or "q4q4"
#
# Returns years x draws, NA where a period that a figure requires is missing.
.annual_figures <- function(values, period, years, multiplicative, target, freq) {

  result <- matrix(NA_real_, length(years), ncol(values))
  for (i in seq_along(years)) {
    current <- years[i] * freq + 0:(freq - 1)
    if (target == "q4q4") {
      current <- current[freq]
    }
    rows <- match(current, period)
    if (anyNA(rows)) {
      next
    }
    level <- colMeans(values[rows, , drop = FALSE])
    if (!multiplicative) {
      result[i, ] <- level
      next
    }
    rows <- match(current - freq, period)
    if (anyNA(rows)) {
      next
    }
    result[i, ] <- 100 * (level / colMeans(values[rows, , drop = FALSE]) - 1)
  }

  return(result)
}

# The annual figures of the data, one column per variable, from the year before
# the first period of the data to the last year the data cover completely. From
# the levels where they are given, which reach as far as they do, and from the
# data of the models otherwise.
.annual_data <- function(original, levels, code, target, scale, freq) {

  source <- if (is.null(levels)) original else levels
  period <- round(stats::time(source) * freq)
  years <- (period[1] %/% freq - 1):((period[length(period)] + 1) %/% freq - 1)

  result <- matrix(NA_real_, length(years), length(code))
  for (j in seq_along(code)) {
    x <- as.numeric(source[, names(code)[j]])
    if (is.null(levels) && code[[j]] >= 4) {
      # The level up to a factor, over the periods the data cover without a gap
      z <- x / scale
      avail <- which(!is.na(z))
      run <- avail[1]:(avail[1] + which(c(diff(avail), 2) != 1)[1] - 1)
      x <- rep(NA_real_, length(z))
      if (code[[j]] == 4) {
        x[run] <- exp(z[run] - z[run[1]])
      } else {
        # The first change is the one out of the period before the run
        x[run] <- exp(cumsum(z[run]))
        if (run[1] > 1) {
          x[run[1] - 1] <- 1
        }
      }
    }
    result[, j] <- .annual_figures(matrix(x, ncol = 1), period, years, code[[j]] >= 4,
                                   target, freq)
  }
  dimnames(result) <- list(NULL, names(code))

  return(stats::ts(result, start = years[1], frequency = 1))
}

# Refuses levels, which do not reproduce the data of the models under the codes
# and the scale, over the periods both cover. A scale of 1 for data in percent,
# or levels of another vintage, would otherwise continue the wrong series.
.check_levels <- function(levels, original, code, scale, freq) {

  for (j in seq_along(code)) {
    var <- names(code)[j]
    z <- .transform_variables_code(as.numeric(levels[, var]), code[[j]])
    if (code[[j]] >= 4) {
      z <- z * scale
    }
    data <- .values_at(original[, var], round(stats::time(levels) * freq), freq)
    both <- !is.na(z) & !is.na(data)
    if (!any(both)) {
      stop("Argument 'levels' does not overlap with the data of variable '", var, "'.")
    }
    if (!isTRUE(all.equal(z[both], data[both], tolerance = 1e-6))) {
      stop("Argument 'levels' does not reproduce the data of variable '", var,
           "' under code ", code[[j]], if (code[[j]] >= 4) paste0(" and scale ", scale),
           ". Check\narguments 'code' and 'scale'.")
    }
  }

  invisible(NULL)
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
