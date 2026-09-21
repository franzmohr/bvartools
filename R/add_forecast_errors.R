#' Add Forecast Errors
#'
#' Generic function used to calculate forecast errors and add them to a model object.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param test_sample a time-series object used as test data. If \code{NULL}
#' (default), the values in \code{data$test$y} of the object are used, which is
#' what a model carries after it has been scored once and what a model read from
#' a file was written with.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_forecast_errors.bvarmodel}},
#' \code{\link{add_forecast_errors.expandingwindow}},
#' \code{\link{add_forecast_errors.modellist}}.
#'
#' @export
add_forecast_errors <- function (object, test_sample = NULL, ...) {
 UseMethod("add_forecast_errors")
}


# The periods of a test sample that a model's forecast covers, one row per
# period and one column per variable, or NULL where the sample does not reach
# them -- which is not an error: a window near the end of a series forecasts
# past what was ever observed, and the windows before it can still be scored.
.align_test_sample <- function(object, test_sample) {

  k <- object[["model"]][["k"]]
  if (k == 1) {
    tsp_test_sample <- stats::tsp(test_sample)
    test_sample <- stats::ts(as.matrix(test_sample), class = c("mts", "ts", "matrix"))
    stats::tsp(test_sample) <- tsp_test_sample
  } else {
    test_sample <- test_sample[, object[["model"]][["endogen"]]]
  }
  test_sample <- stats::na.omit(test_sample)

  # Determine when the forecasts start
  tsp_train <- stats::tsp(object[["data"]][["train"]][["y"]])

  # A period of the test sample at another frequency would be matched by its
  # time alone: 2020 is a year and the first quarter of 2020 alike
  if (!isTRUE(all.equal(stats::tsp(test_sample)[3], tsp_train[3]))) {
    stop("Argument 'test_sample' has ", frequency_name(stats::tsp(test_sample)[3]),
         " data, but the forecasts of the model are ", frequency_name(tsp_train[3]), ".")
  }
  forecast_starts_at <- tsp_train[2] + 1 / tsp_train[3]

  if (!(forecast_starts_at %in% stats::time(test_sample))) {
    return(NULL)
  }

  test_sample <- stats::window(test_sample, start = forecast_starts_at)
  h <- min(object[["model"]][["h"]], nrow(test_sample))

  return(as.matrix(test_sample[1:h, , drop = FALSE]))
}


# The realised values a model carries in data$test$y, checked against what it
# forecast. Already aligned: they are the periods of the horizon and nothing
# else, which is what is written to /data/test/y of a model file.
.realised_values <- function(object) {

  y <- object[["data"]][["test"]][["y"]]
  # Forecasts aggregated to annual figures carry the realised values from the
  # start, and a window without them forecasts years the data do not cover yet
  if (is.null(y) && !is.null(object[["model"]][["aggregation"]])) {
    return(NULL)
  }
  if (is.null(y)) {
    stop("No test sample was given and the object carries none in data$test$y. ",
         "Pass 'test_sample', or score a model that was written to a file after ",
         "'add_forecast_errors' had put the values it was scored against into it.")
  }

  y <- as.matrix(y)
  k <- object[["model"]][["k"]]
  h <- object[["model"]][["h"]]
  if (ncol(y) != k) {
    stop("Element data$test$y has ", ncol(y), " columns, but the model has ", k,
         " endogenous variables. It is one row per period and one column per variable.")
  }
  if (nrow(y) > h) {
    stop("Element data$test$y holds ", nrow(y), " periods, more than the ", h,
         " this model forecasts.")
  }

  return(y)
}
