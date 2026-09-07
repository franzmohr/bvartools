#' Model Selection Criteria
#'
#' Calculates out-of-sample statistics for an object of class 'externalforecast'.
#'
#' @param object an object of class 'externalforecast'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class 'selcrit'.
#'
#' @details
#' External forecasts are not estimated, so only out-of-sample statistics are
#' calculated. Since they are point forecasts, each publication contributes a single
#' value to the statistics of a variable and forecast horizon.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' # Create model
#' model <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "none",
#'                           error = "gamma", iterations = 10, burnin = 2)
#' # Chosen number of iterations and burn-in draws should be much higher.
#'
#' model <- use_expanding_window(model, start = 2007)
#'
#' # Artificial external forecasts
#' fcst <- expand.grid(origin = c(2007, 2007.25),
#'                     h = 1:2,
#'                     variable = c("Dp", "r"),
#'                     stringsAsFactors = FALSE)
#' fcst[["period"]] <- fcst[["origin"]] + fcst[["h"]] / 4
#' fcst[["value"]] <- 0
#'
#' ext <- create_external_forecast(fcst, model, n_ahead = 4)
#'
#' # Calculate forecast errors
#' ext <- add_forecast_errors(ext, test_sample = us_macrodata)
#'
#' # Calculate selection criteria
#' selection_criteria(ext)
#'
#' @export
#' @method selection_criteria externalforecast
selection_criteria.externalforecast <- function(object, ci = 0.95, ...) {

  if (ci < 0 | ci > 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low

  k <- object[[1]][["model"]][["k"]]
  h <- object[[1]][["model"]][["h"]]
  varnames <- object[[1]][["model"]][["endogen"]]
  max_n_columns <- k * h

  fcst_errors <- lapply(object, get_forecast_errors, ...)

  # The forecast errors of a publication can cover fewer periods than the maximum
  # forecast horizon, if the test data end before the last forecasted period
  errors <- NULL
  for (i in 1:length(fcst_errors)) {
    if (!is.null(fcst_errors[[i]])) {
      if (ncol(fcst_errors[[i]]) < max_n_columns) {
        empty_matrix <- matrix(NA_real_, nrow(fcst_errors[[i]]),
                               max_n_columns - ncol(fcst_errors[[i]]))
        fcst_errors[[i]] <- cbind(fcst_errors[[i]], empty_matrix)
      }
      errors <- rbind(errors, fcst_errors[[i]])
    }
  }

  if (is.null(errors)) {
    stop("Object does not contain any forecast errors. You might want to use function\n'add_forecast_errors' before this function.")
  }

  result <- NULL
  result[["model"]] <- object[[1]][["model"]]

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
  # The classes of the publications are maintained, so that methods, which use the
  # model specifications, can be dispatched on them
  class(result) <- c("selcrit", class(object[[1]]))

  return(result)
}
