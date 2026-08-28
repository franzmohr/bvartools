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
#' @return A list of class 'expandingwindow'.
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
    
    result[[i]] <- temp
    rm(temp)
  }
  
  class(result) <- append("expandingwindow", class(result))
  
  return(result)
}