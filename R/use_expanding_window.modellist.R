#' Expanding Window Estimation
#' 
#' Creates objects for expanding window posterior simulation.
#' 
#' @param object a list of objects containing model specifications and input data.
#' Usually, the output of a call to a model creation function such as
#' \code{\link{create_bvarmodel}} or \code{\link{create_bvecmodel}}.
#' @param start the start period of the prediction of the first iteration of the
#' expanding window approach.
#' @param ... arguments passed forward to method.
#' 
#' @returns A list of class 'modellist', which consists of objects of class
#' 'expandingwindow'.
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
#'                           p = 1:4,
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
#' @method use_expanding_window modellist
use_expanding_window.modellist <- function(object, start, ...) {
  
  result <- list()
  for (i in 1:length(object)) {
    result[[i]] <- use_expanding_window(object[[i]], start = start, ...)
  }
  
  class(result) <- append("modellist", class(result))
  
  return(result)
}