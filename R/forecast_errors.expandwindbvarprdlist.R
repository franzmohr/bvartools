#' Forecast Errors
#'
#' Calculates forecast errors for a list of Bayesian models.
#'
#' @param object an object of class 'expandwindbvarprdlist', usually, the
#' result of a call to \code{\link{predict.expandwindbvarprdlist}}.
#' @param test_sample a time-series object used as test data.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'fcsterrorlist'.
#' 
#' @examples
#' 
#' data("us_macrodata")
#' 
#' # AR(1) models as benchmark
#' model <- create_bvarmodel(data = us_macrodata,
#'                           p = 1,
#'                           deterministic = "none",
#'                           error = "gamma",
#'                           iterations = 10,
#'                           burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#' 
#' # Obtain objects for expanding window estimation
#' model <- use_expanding_window(model, start = 2007)
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(shape = 3, rate = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' post <- draw_posterior(model)
#' 
#' # Predict
#' prd <- predict(post, n_ahead = 4)
#' 
#' # Forecast errors
#' fe <- forecast_errors(prd, test_sample = us_macrodata)
#' 
#' 
#' @export
forecast_errors.expandwindbvarprdlist <- function(object, test_sample, ...){
  
  object <- lapply(object, forecast_errors, test_sample = test_sample)
  
  # Get n_ahead
  n_ahead <- max(unlist(lapply(object, function(x) {dim(x)[1]})))
  
  # Number variables
  k <- unique(unlist(lapply(object, function(x) {dim(x)[2]})))
  
  # Number of draws
  draws <- unique(unlist(lapply(object, function(x) {dim(x)[3]})))
  
  # Get time
  tsp <- do.call("rbind", lapply(object, function(x) {attr(x, "tsp")}))
  tsp <- unique(c(tsp[, 1:2]))
  tsp <- tsp[order(tsp)]
  
  result <- list()
  for (i in 1:k) {
    result[[i]] <- array(NA_real_,
                         dim = c(length(tsp), n_ahead, draws * length(object)),
                         dimnames = list(as.character(tsp), NULL, NULL))
    class(result[[i]]) <- append("fcsterror", class(result[[i]]))
  }
  vars <- dimnames(object[[1]])[[2]]
  names(result) <- vars
  
  for (var in vars) {
    for (i in 1:length(object)) {
      if (is.null(object[[i]])) {
        next
      }
      pos_time <- dimnames(object[[i]])[[1]]
      dims_i <- dim(object[[i]])
      for (h in 1:n_ahead) {
        if (h <= dims_i[1]) {
          result[[var]][pos_time[h], h, (i - 1) * draws + 1:draws] <- object[[i]][h, var,] 
        }
      }
    }
  }
  
  class(result) <- append("bvarfcsterror", class(result))
  
  return(result)
}