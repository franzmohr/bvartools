#' Expanding Window Estimation
#' 
#' Creates objects for expanding window posterior simulation.
#' 
#' @param object an object of class 'bvecmodel' containing model specifications
#' and input data. Usually, the output of a call to  \code{\link{create_bvecmodel}}.
#' @param start the start period of the prediction of the first iteration of the
#' expanding window approach.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'expandingwindow'.
#' 
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Create multiple model objects for expanding window
#' model <- use_expanding_window(model, start = 1998)
#' 
#' @export
#' @method use_expanding_window bvecmodel
use_expanding_window.bvecmodel <- function(object, start, ...) {

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
    
    temp[["data"]][["train"]][["y"]] <- stats::window(temp[["data"]][["train"]][["y"]], end = time_y[pos_end[i]])
    if (!is.null(temp[["data"]][["train"]][["w"]])) {
      temp[["data"]][["train"]][["w"]] <- stats::window(temp[["data"]][["train"]][["w"]], end = time_y[pos_end[i]])
    }
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