#' Add Forecast Input Data
#'
#' Generates and adds data matrices for forecast simulation to the elements of
#' an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow' containing objects that can
#' be passed on to their respective \code{add_forecast_input} method.
#' @param ... arguments passed forward to method. Argument \code{y_start} of
#' \code{\link{add_forecast_input.bvarmodel}} is refused, since every window forecasts from
#' the end of its own estimation sample.
#'
#' @return The object in \code{object} with forecast input added to each of its models, as
#' described in \code{\link{add_forecast_input.bvarmodel}}.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' orig <- diff(log(e1)) * 100
#' train <- window(orig, end = c(1982, 2))
#' 
#' 
#' # Create model
#' model <- create_bvarmodel(data = train, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' model <- use_expanding_window(model, start = 1982.25)
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#' 
#' @export
#' @method add_forecast_input expandingwindow
add_forecast_input.expandingwindow <- function(object, ...){

  # Every window forecasts the periods after its own sample, which is what the
  # exercise compares. One set of start values would make them all start in the
  # same period.
  if (!is.null(list(...)[["y_start"]])) {
    stop("Argument 'y_start' is not available for an 'expandingwindow' object: each window ",
         "forecasts from the end of its own estimation sample. Give it to add_forecast_input() ",
         "of a single model instead, such as the last window, object[[length(object)]].",
         call. = FALSE)
  }

  for (i in 1:length(object)) {
    object[[i]] <- add_forecast_input(object[[i]], ...)
  }
  
  return(object)
}
