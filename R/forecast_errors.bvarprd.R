#' Forecast Errors
#'
#' Calculates forecast errors for models of class 'bvarprd'.
#'
#' @param object an object of class 'bvarprd', usually, the
#' result of a call to \code{\link{predict.bvarmodel}}.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'fcsterror'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1_long <- e1
#' e1 <- window(e1_long, end = c(1976, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 0, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Generate forecasts
#' bvar_pred <- predict(object, n_ahead = 10)
#' 
#' # Forecast errors
#' fe <- forecast_errors(bvar_pred, test_sample = e1_long)
#' 
#' 
#' @export
forecast_errors.bvarprd <- function(object, test_sample, ...){
  
  dims <- dim(object[["fcst"]])
  n_ahead <- dims[1]
  k <- dims[2]
  draws <- dims[3]
  tsp_fcst <- attr(object[["fcst"]], "tsp")
  
  if (any(stats::time(test_sample) >= tsp_fcst[1])) {
    pos_orig <- which(stats::time(test_sample) >= tsp_fcst[1])
    n_ahead <- min(length(pos_orig), n_ahead)
    pos_orig <- pos_orig[1:n_ahead]
    pos <- 1:n_ahead
    test_sample <- matrix(test_sample[pos_orig,], n_ahead)
    
    result <- object[["fcst"]] * NA_real_
    for (draw in 1:draws) {
      fc_draw <- stats::ts(matrix(object[["fcst"]][pos,,draw], n_ahead), start = tsp_fcst[1], frequency = tsp_fcst[3])
      result[pos,,draw] <- test_sample - fc_draw
    } 
    
    class(result) <- c("fcsterror", "array")
  } else {
    result <- NULL
  }
  
  return(result)
}