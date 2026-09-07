#' @include create_external_forecast.R
#'
#' @param x an object of class 'externalforecast'.
#' @param thin an integer specifying the thinning interval between successive draws.
#' @param digits the number of significant digits.
#' @param ... arguments passed forward to method.
#'
#' @details
#' External forecasts do not have to be estimated, so the functions, which add priors,
#' initial values or posterior draws to a model are without effect for objects of
#' class 'externalforecast'. This allows to combine them with model objects in a list
#' of class 'modellist' and to apply the usual workflow to all elements of that list.
#'
#' @rdname create_external_forecast
#' @export
add_priors.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
add_initial_values.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
add_posterior_coefficients.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
add_posterior_loglik.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
add_forecast_input.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
add_posterior_forecasts.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
minnesota_prior.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
ssvs_prior.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
inclusion_prior.externalforecast <- function(object, ...) {
  return(object)
}

#' @rdname create_external_forecast
#' @export
thin.externalforecast <- function(x, thin = 10, ...) {
  return(x)
}

#' @rdname create_external_forecast
#' @export
print.externalforecast <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("External forecasts\n\n")

  cat("Variables: ", paste0(x[[1]][["model"]][["endogen"]], collapse = ", "), "\n", sep = "")
  cat("Forecast horizon: ", x[[1]][["model"]][["h"]], "\n", sep = "")
  cat("Publications: ", length(x), "\n\n", sep = "")

  # The periods are formatted as characters, because 'print' would otherwise round
  # them to the number of significant digits
  result <- data.frame("Publication" = format(unlist(lapply(x, function(y) {
                         y[["model"]][["origin"]]
                       })), trim = TRUE),
                       "Training sample ends" = format(unlist(lapply(x, function(y) {
                         stats::tsp(y[["data"]][["train"]][["y"]])[2]
                       })), trim = TRUE),
                       check.names = FALSE, stringsAsFactors = FALSE)
  print(result, digits = digits, row.names = FALSE, ...)

  return(invisible(x))
}
