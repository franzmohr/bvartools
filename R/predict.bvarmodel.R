#' Predict Method for Objects of Class bvar
#' 
#' Forecasting a Bayesian VAR object of class 'bvarmodel'.
#' 
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with deterministic data. If not
#' specified, the function will try to identify the deterministic terms
#' automatically. If this is not successful, an error message we be returned.
#' @param exogen a time-series object with unmodeled, non-deterministic data.
#' See 'Details'.
#' @param ... additional arguments.
#' 
#' @details For the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + \sum_{i = 0}^{s} B_{i} x_{t-i} + C D_t + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)} the function produces \code{n_ahead} forecasts.
#' 
#' Data provided in argument \code{exogen} will be prepared according to the
#' model specifications in argument \code{object} and joined with the prediction
#' data set using the time stamp as a key. Therefore, \code{exogen} should
#' only contain the series that should be used. It is \emph{not} necessary to add further
#' columns with lags of the time series.
#' 
#' @return An array of class 'bvarprd'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
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
#' model <- add_posterior_coefficients(model)
#' 
#' # Calculate forecasts
#' pred <- predict(model)
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
predict.bvarmodel <- function(object, n_ahead = 10, ...) {
  
  # Prepare input
  if (is.null(object[["model"]][["h"]])) {
    stop("Missing specification of h in object$model$h. You might want to use\nfunction 'add_forecast_input' and then 'add_posterior_forecasts'\nbefore using this function.")
  }
  
  if (is.null(object[["posterior"]][["forecast"]])) {
    stop("Missing element object$posterior$forecast. You might want to use\n'add_posterior_forecasts'\nbefore this function.")
  }
  
  if (object[["model"]][["h"]] < n_ahead) {
    warning("Argument 'n_ahead' is larger than the value in object$model$h.\nLimiting the output to the latter.")
    n_ahead <- object[["model"]][["h"]]
  }
  
  tt <- nrow(object[["data"]][["train"]][["y"]])
  varnames <- object[["model"]][["endogen"]]
  draws <- nrow(object[["posterior"]][["forecast"]])
  
  tsp_temp <- stats::tsp(object[["data"]][["train"]][["y"]])
  prd_time <- tsp_temp[2] + 1:object[["model"]][["h"]] / tsp_temp[3]
  
  result <- array(NA_real_, dim = c(object[["model"]][["h"]], length(varnames), draws))
  dimnames(result) <- list(as.character(prd_time), varnames, NULL)
  attr(result, "tsp") <- c(min(prd_time), max(prd_time), tsp_temp[3])
  
  for (i in 1:draws) {
    result[,,i] <- t(matrix(object[["posterior"]][["forecast"]][i,], object[["model"]][["k"]]))
  }
  
  result <- list(fcst = result,
                 y = object[["data"]][["train"]][["y"]])
  
  class(result) <- c("bvarprd", "list")
  return(result)
}
