#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'selcrit'.
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
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Add data used for forecast calculation
#' model <- add_forecast_input(model, n_ahead = 4)
#' 
#' # Add forecasts
#' model <- add_posterior_forecasts(model)
#' 
#' # Add forecast errors
#' model <- add_forecast_errors(model, test_sample = orig)
#' 
#' # Calculate selection criteria
#' sel <- selection_criteria(model)
#' sel
#' 
#' 
#' @export
selection_criteria.expandingwindow <- function(object, ci = 0.95, ...){
  
  if (ci < 0 | ci > 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  
  
  # Model information
  k <- object[[1]][["model"]][["k"]]
  h <- object[[1]][["model"]][["h"]]
  varnames <- object[[1]][["model"]][["endogen"]]
  max_n_columns <- k * h
  
  fcst_errors <- lapply(object, get_forecast_errors, ...)
  
  errors <- NULL
  for (i in 1:length(fcst_errors)) {
    if (!is.null(fcst_errors[[i]])) {
      if (ncol(fcst_errors[[i]]) == max_n_columns) {
        errors <- rbind(errors, fcst_errors[[i]])
      } else {
        emtpy_matrix <- matrix(NA_real_, nrow(fcst_errors[[i]]), max_n_columns - ncol(fcst_errors[[i]]))
        errors <- rbind(errors, cbind(fcst_errors[[i]], emtpy_matrix))
      }
    }
  }
  
  
  result <- NULL
  result[["model"]] <- object[[1]][["model"]]
  
  # In-sample
  in_sample <- selection_criteria(object[[length(object)]], ci = ci)
  for (i in c("LL", "AIC", "BIC", "HQ")) {
    result[[i]] <- in_sample[[i]]
  }
  
  
  # Out-ofsample
  
  # Forecast errors
  result[["FE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
  names(result[["FE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
  result[["FE"]][, "variable"] <- rep(varnames, h)
  result[["FE"]][, "h"] <- rep(1:h, each = k)
  result[["FE"]][, "mean"] <- apply(errors, 2, mean, na.rm = TRUE)
  result[["FE"]][, "median"] <- apply(errors, 2, stats::median, na.rm = TRUE)
  result[["FE"]][, "qlower"] <- apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE)
  result[["FE"]][, "qupper"] <- apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)
  
  # Absolute errors
  errors <- abs(errors)
  result[["AFE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
  names(result[["AFE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
  result[["AFE"]][, "variable"] <- rep(varnames, h)
  result[["AFE"]][, "h"] <- rep(1:h, each = k)
  result[["AFE"]][, "mean"] <- apply(errors, 2, mean, na.rm = TRUE)
  result[["AFE"]][, "median"] <- apply(errors, 2, stats::median, na.rm = TRUE)
  result[["AFE"]][, "qlower"] <- apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE)
  result[["AFE"]][, "qupper"] <- apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)
  
  # Squared errors
  errors <- errors^2
  result[["RSFE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
  names(result[["RSFE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
  result[["RSFE"]][, "variable"] <- rep(varnames, h)
  result[["RSFE"]][, "h"] <- rep(1:h, each = k)
  result[["RSFE"]][, "mean"] <- sqrt(apply(errors, 2, mean, na.rm = TRUE))
  result[["RSFE"]][, "median"] <- sqrt(apply(errors, 2, stats::median, na.rm = TRUE))
  result[["RSFE"]][, "qlower"] <- sqrt(apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE))
  result[["RSFE"]][, "qupper"] <- sqrt(apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)) 
  
  attr(result, "ci") <- c(paste0(ci_low * 100, "%"), paste0(ci_high * 100, "%"))
  # The classes of the models of the estimation windows are maintained, so that
  # methods, which use the model specifications, can be dispatched on them
  class(result) <- c("selcrit", class(object[[1]]))
  
  return(result)
}