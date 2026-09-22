#' Expanding Window Estimation
#' 
#' Creates objects for expanding window posterior simulation.
#' 
#' @param object an object of class 'bvarmodel' containing model specifications
#' and input data. Usually, the output of a call to  \code{\link{create_bvarmodel}}.
#' @param start the start period of the prediction of the first iteration of the
#' expanding window approach.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'expandingwindow' with one object of class 'bvarmodel' per
#' window. The training sample of the first ends in the period before \code{start}, and
#' each further window adds one period. Posterior draws, starting values and forecast input
#' that \code{object} already carries belong to the whole sample and are not copied into
#' the windows; a warning says so.
#'
#' @examples
#' 
#' data("us_macrodata")
#' 
#' # Starting period of the forecasting exercise
#' start_period <- 2007
#' 
#' # Create model
#' model <- create_bvarmodel(data = us_macrodata,
#'                           p = 1,
#'                           deterministic = "none",
#'                           seasonal = FALSE,
#'                           tvp = FALSE,
#'                           error = "gamma",
#'                           iterations = 10,
#'                           burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Create multiple model objects for expanding window
#' model <- use_expanding_window(model, start = start_period)
#' 
#' @family model set-up
#' @export
#' @method use_expanding_window bvarmodel
use_expanding_window.bvarmodel <- function(object, start, ...) {

  k <- object[["model"]][["k"]]
  y <- object[["data"]][["train"]][["y"]]
  test <- stats::window(y, start = start)
  time_y <- stats::time(y)
  time_train <- time_y[time_y < min(stats::time(test))]
  nobs_train_min <- length(time_train)
  nobs_train_max <- length(time_y)
  pos_end <- nobs_train_min:nobs_train_max
  
  object <- .clear_for_windows(object)

  # Produce individual models with incrementally increasing estimation horizons
  result <- list()
  for (i in 1:length(pos_end)) {

    temp <- object
    
    # Trim data
    temp[["data"]][["train"]][["y"]] <- stats::window(temp[["data"]][["train"]][["y"]], end = time_y[pos_end[i]])
    if (!is.null(temp[["data"]][["train"]][["x"]])) {
      temp[["data"]][["train"]][["x"]] <- stats::window(temp[["data"]][["train"]][["x"]], end = time_y[pos_end[i]])
    }
    if (!is.null(temp[["data"]][["train"]][["z"]])) {
      temp[["data"]][["train"]][["z"]] <- temp[["data"]][["train"]][["z"]][1:(k * pos_end[i]), ]
    }
    
    # Every window is simulated on its own, so a seed the model already has is
    # counted up from window to window rather than shared.
    if (!is.null(temp[["model"]][["seed"]])) {
      temp[["model"]][["seed"]] <- .offset_seed(temp[["model"]][["seed"]], i - 1)
    }

    result[[i]] <- temp
    rm(temp)
  }
  
  class(result) <- append("expandingwindow", class(result))
  
  return(result)
}


# What a model carries that was estimated on, or built for, the whole sample,
# removed before the model is copied into its windows. The posterior is the one
# that matters most: a window that kept it would forecast and be scored with
# draws that had seen the periods it is evaluated on, and nothing downstream
# could tell. The starting values can be paths as long as the whole sample, and
# the forecast input continues the whole sample rather than a window. Priors
# are kept -- setting them before splitting the sample is how they are meant to
# be given once for every window.
.clear_for_windows <- function(object) {
  dropped <- character(0)
  if (!is.null(object[["posterior"]])) {
    object[["posterior"]] <- NULL
    dropped <- c(dropped, "the posterior draws")
  }
  if (!is.null(object[["initial"]])) {
    object[["initial"]] <- NULL
    dropped <- c(dropped, "the starting values")
  }
  if (!is.null(object[["data"]][["forecast"]]) || !is.null(object[["data"]][["test"]])) {
    object[["data"]][["forecast"]] <- NULL
    object[["data"]][["test"]] <- NULL
    object[["model"]][["h"]] <- NULL
    dropped <- c(dropped, "the forecast input and test sample")
  }
  if (length(dropped) > 0) {
    warning("Argument 'object' already carries ", paste(dropped, collapse = ", "),
            ", which belong to the whole sample and were not copied into the windows. ",
            "Add them to the result instead, for example with add_initial_values(), ",
            "add_posterior_coefficients() and add_forecast_input().", call. = FALSE)
  }
  object
}
