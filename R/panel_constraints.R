#' Weights of a Temporal Aggregate
#'
#' The weights with which a series observed at a lower frequency aggregates the
#' high-frequency periods of a model, for argument \code{aggregate} of
#' \code{\link{create_bvarmodel}}.
#'
#' @param type a character naming the aggregation: \code{"average"} for a
#' series that is the mean of the periods it covers, such as a quarterly
#' unemployment rate in a monthly model, \code{"sum"} for a flow that is their
#' total, and \code{"growth"} for the growth rate of an average, following
#' Mariano and Murasawa (2003). See 'Details'.
#' @param n an integer, the number of high-frequency periods one low-frequency
#' period covers: 3 for quarters of months, 12 for years of months, 4 for years
#' of quarters.
#'
#' @details An observation of the low-frequency series in period \eqn{t} is
#' \deqn{x_t = \sum_{j = 1}^{L} w_j y_{t - L + j},}
#' where \eqn{y_t} is the high-frequency series the model contains and
#' \eqn{w_1, \dots, w_L} are the weights returned, oldest first.
#'
#' For \code{"average"} they are \eqn{n} weights of \eqn{1/n}, and for
#' \code{"sum"} \eqn{n} weights of one. For \code{"growth"} the high-frequency
#' series is the growth rate -- the log difference -- of a series whose average
#' the low-frequency series is the growth rate of. Approximating the log of an
#' average by the average of the logs, that growth rate is a triangular
#' combination of \eqn{2n - 1} high-frequency growth rates, which for quarters
#' of months is \eqn{(1, 2, 3, 2, 1) / 3}.
#'
#' @return A numeric vector of weights, oldest period first.
#'
#' @examples
#' aggregation_weights("growth", 3)
#'
#' @references
#'
#' Mariano, R. S., & Murasawa, Y. (2003). A new coincident index of business
#' cycles based on monthly and quarterly series. \emph{Journal of Applied
#' Econometrics, 18}(4), 427--443. \doi{10.1002/jae.695}
#'
#' @family model set-up
#' @export
aggregation_weights <- function(type = c("average", "sum", "growth"), n = 3) {

  type <- match.arg(type)
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n %% 1 != 0) {
    stop("Argument 'n' must be a whole number of at least one.")
  }

  switch(type,
         "average" = rep(1 / n, n),
         "sum" = rep(1, n),
         "growth" = c(seq_len(n), rev(seq_len(n - 1))) / n)
}


# The panel of a model estimated from what was observed of it rather than from
# the periods where every series was. See create_bvarmodel().
#
# Checks 'aggregate' and 'soft' against the series in 'data' and returns them
# with 'aggregate' as a named list of weights for the series it names.
.check_panel_arguments <- function(data, aggregate, soft) {

  series <- dimnames(data)[[2]]

  if (!is.null(aggregate)) {
    if (!is.list(aggregate) || is.null(names(aggregate)) || any(names(aggregate) == "")) {
      stop("Argument 'aggregate' must be a named list of weights, one element per aggregated series.")
    }
    unknown <- setdiff(names(aggregate), series)
    if (length(unknown) > 0) {
      stop("Argument 'aggregate' names series that are not in 'data': ",
           paste(unknown, collapse = ", "), ".")
    }
    for (i in names(aggregate)) {
      w <- aggregate[[i]]
      if (!is.numeric(w) || length(w) == 0 || anyNA(w) || any(!is.finite(w)) || sum(w) == 0) {
        stop("The weights of series '", i, "' in argument 'aggregate' must be finite numbers ",
             "that do not sum to zero.")
      }
    }
  }

  if (!is.null(soft)) {
    if (!is.character(soft) || anyNA(soft) || anyDuplicated(soft)) {
      stop("Argument 'soft' must be a character vector of distinct series names.")
    }
    unknown <- setdiff(soft, series)
    if (length(unknown) > 0) {
      stop("Argument 'soft' names series that are not in 'data': ",
           paste(unknown, collapse = ", "), ".")
    }
  }

  list("aggregate" = aggregate, "soft" = soft)
}

# 'data' with every gap filled, so that the lags and the regressors can be built
# as for a complete panel. What stands in a gap is a starting value and nothing
# more: the sampler redraws every period that was not observed. It is the linear
# interpolation of what was, carried flat past either end, and for an aggregated
# series the interpolation of its observations divided by the sum of their
# weights, which puts the start on the scale of the series the model contains.
.fill_panel <- function(data, aggregate) {

  series <- dimnames(data)[[2]]
  filled <- data
  index <- seq_len(nrow(data))
  for (j in seq_along(series)) {
    values <- as.numeric(data[, j])
    observed <- which(!is.na(values))
    if (length(observed) == 0) {
      stop("Series '", series[j], "' is never observed.")
    }
    if (!is.null(aggregate[[series[j]]])) {
      values <- values / sum(aggregate[[series[j]]])
    }
    if (length(observed) == 1) {
      filled[, j] <- values[observed]
    } else {
      filled[, j] <- stats::approx(observed, values[observed], xout = index, rule = 2)$y
    }
  }

  filled
}

# What was observed of the estimation sample 'y', as the constraint set
# data$train$constraints holds: one row per observation, the six elements of
# /data/train/constraints of a model file, positions counted from one and hard
# rows in group 0. NULL when the sample is observed whole.
#
# Observations are placed by date rather than by position, so that the lags --
# and exogenous series that start later -- shortening the sample leave them
# where they belong. An observation of an aggregated series reaches back over
# as many periods as it has weights; one that reaches before the sample is
# dropped, since the periods it would constrain are not estimated.
.panel_constraints <- function(original, y, aggregate, soft) {

  series <- dimnames(original)[[2]]
  freq <- stats::frequency(y)
  start <- stats::tsp(y)[1]
  tt <- nrow(y)
  periods <- round((as.numeric(stats::time(original)) - start) * freq) + 1

  if (is.null(aggregate) && !anyNA(original[periods >= 1 & periods <= tt, ])) {
    return(NULL)
  }

  value <- numeric(0)
  group <- numeric(0)
  row <- numeric(0)
  period <- numeric(0)
  variable <- numeric(0)
  weight <- numeric(0)
  dropped <- 0L
  r <- 0

  for (j in seq_along(series)) {
    w <- aggregate[[series[j]]]
    if (is.null(w)) {
      w <- 1
    }
    lag <- length(w) - 1
    g <- if (series[j] %in% soft) match(series[j], soft) else 0
    for (i in which(!is.na(original[, j]))) {
      last <- periods[i]
      if (last > tt || last < 1) {
        next
      }
      if (last - lag < 1) {
        dropped <- dropped + 1L
        next
      }
      r <- r + 1
      value <- c(value, unname(as.numeric(original[i, j])))
      group <- c(group, g)
      row <- c(row, rep(r, length(w)))
      period <- c(period, (last - lag):last)
      variable <- c(variable, rep(j, length(w)))
      weight <- c(weight, w)
    }
  }

  if (dropped > 0) {
    message(dropped, " observation(s) of an aggregated series reach before the estimation ",
            "sample and are left out.")
  }

  # Groups are numbered from one with none left out, so a series named in
  # 'soft' that contributes no row renumbers the ones after it.
  used <- sort(unique(group[group > 0]))
  group[group > 0] <- match(group[group > 0], used)

  list("value" = value, "group" = group, "row" = row,
       "period" = period, "variable" = variable, "weight" = weight)
}
