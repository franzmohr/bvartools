#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class 'expandingwindow'.
#'
#' @param object an object of class 'expandingwindow'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The in-sample criteria are those of the last window, the one estimated on
#' the most data, and the forecast error statistics those of all windows.
#'
#' If the windows hold the draws of \code{\link{add_predictive_loglik}},
#' criterion \code{"LPL"} is the log predictive likelihood: the sum over the
#' windows of the log of the mean of the draws of the predictive density of the
#' observation that the next window adds. Its band is the normal interval of
#' probability \code{ci} around the sum with standard error
#' \eqn{\sqrt{n \mathrm{Var}(lpd_t)}} over the \eqn{n} periods, as for LOOIC.
#' Attribute \code{"terms"} holds the log predictive density of each period and
#' its numerical standard error, attribute \code{"nse"} the numerical standard
#' error of the sum.
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

  # Each half only if its draws exist. An expanding window is estimated for its
  # log-likelihood as often as for its forecasts, and the out-of-sample half
  # used to be built unconditionally, from a NULL matrix of errors, which
  # stopped the whole call on a window that had no forecasts. The in-sample
  # half is that of the last window, the one estimated on the most data.
  last_window <- object[[length(object)]]
  use_ll <- !is.null(last_window[["posterior"]][["loglik"]])
  use_fe <- !is.null(errors)
  predictive <- Filter(Negate(is.null), lapply(object, function(x) x[["predictive"]]))
  use_lpl <- length(predictive) > 0
  if (!use_ll & !use_fe & !use_lpl) {
    stop("Model object must contain at least either posterior draws of the log-likelihood, ",
         "predictive log-likelihoods or forecast errors.")
  }

  # Log predictive likelihood
  if (use_lpl) {
    result[["LPL"]] <- .lpl_entry(
      lapply(predictive, function(p) p[["loglik"]]),
      vapply(predictive, function(p) as.numeric(p[["period"]]), numeric(1)),
      ci_low, ci_high)
  }

  # In-sample
  if (use_ll) {
    in_sample <- selection_criteria(last_window, ci = ci)
    for (i in c("LL", "AIC", "BIC", "HQ", "WAIC", "LOOIC")) {
      result[[i]] <- in_sample[[i]]
    }
  }

  # Out-of-sample
  if (use_fe) {

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
  }

  attr(result, "ci") <- c(paste0(ci_low * 100, "%"), paste0(ci_high * 100, "%"))
  # The classes of the models of the estimation windows are maintained, so that
  # methods, which use the model specifications, can be dispatched on them
  class(result) <- c("selcrit", class(object[[1]]))
  
  return(result)
}