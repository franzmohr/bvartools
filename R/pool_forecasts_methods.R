#' @include pool_forecasts.R
#'
#' @param object,x an object of class \code{"forecastpool"} or \code{"poolwindow"}.
#' @param thin an integer specifying the thinning interval between successive draws.
#' @param seed not used, since a pool is not simulated.
#' @param digits the number of significant digits.
#'
#' @details
#' A pool is formed from draws its members already hold, so the functions that
#' add priors, initial values, posterior draws, forecasts or predictive
#' log-likelihoods to a model return it unchanged. This allows to combine a pool
#' with its members in a list of class \code{"modellist"} and to apply the usual
#' workflow to all elements of that list.
#'
#' @rdname pool_forecasts
#' @export
add_priors.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_initial_values.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_seed.forecastpool <- function(object, seed, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_coefficients.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_loglik.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_forecast_input.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_forecasts.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_predictive_loglik.forecastpool <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
thin.forecastpool <- function(x, thin = 10, ...) {
  return(x)
}

#' @rdname pool_forecasts
#' @export
add_priors.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_initial_values.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_seed.poolwindow <- function(object, seed, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_coefficients.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_loglik.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_forecast_input.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_posterior_forecasts.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
add_predictive_loglik.poolwindow <- function(object, ...) {
  return(object)
}

#' @rdname pool_forecasts
#' @export
thin.poolwindow <- function(x, thin = 10, ...) {
  return(x)
}

#' @rdname pool_forecasts
#' @export
write_to_hdf5.forecastpool <- function(object, ...) {
  stop("A pool is formed from the draws of its members and is not written to a file. ",
       "Write the members, and form the pool again with 'pool_forecasts' after reading them back.")
}

#' @rdname pool_forecasts
#' @export
write_to_hdf5.poolwindow <- function(object, ...) {
  stop("A pool is formed from the draws of its members and is not written to a file. ",
       "Write the members, and form the pool again with 'pool_forecasts' after reading them back.")
}

#' @rdname pool_forecasts
#' @export
print.forecastpool <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  first <- x[[1]][["model"]]
  cat("Pooled forecasts, equal weights\n\n")
  cat("Members: ", paste0(first[["members"]], collapse = ", "), "\n", sep = "")
  cat("Variables: ", paste0(first[["endogen"]], collapse = ", "), "\n", sep = "")
  cat("Forecast horizon: ", first[["h"]], "\n", sep = "")
  cat("Draws per member: ", first[["draws"]], "\n", sep = "")
  ends <- vapply(x, .pool_sample_end, numeric(1))
  freq <- stats::frequency(x[[1]][["data"]][["train"]][["y"]])
  cat("Forecasts: ", length(x), ", from estimation samples ending in ",
      .pool_period_label(min(ends), freq), " to ", .pool_period_label(max(ends), freq),
      "\n", sep = "")

  return(invisible(x))
}

#' @rdname pool_forecasts
#' @export
print.poolwindow <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("Pooled forecast, equal weights\n\n")
  cat("Members: ", paste0(x[["model"]][["members"]], collapse = ", "), "\n", sep = "")
  cat("Variables: ", paste0(x[["model"]][["endogen"]], collapse = ", "), "\n", sep = "")
  cat("Forecast horizon: ", x[["model"]][["h"]], "\n", sep = "")
  cat("Draws per member: ", x[["model"]][["draws"]], "\n", sep = "")
  cat("Estimation sample ends in ",
      .pool_period_label(.pool_sample_end(x), stats::frequency(x[["data"]][["train"]][["y"]])),
      "\n", sep = "")

  return(invisible(x))
}

# A time of a ts object as a period label, e.g. 2017.75 as "2017Q4"
.pool_period_label <- function(time, freq) {
  if (freq == 1) {
    return(as.character(round(time)))
  }
  .format_period(round(time * freq), freq)
}
