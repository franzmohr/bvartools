#' Posterior Simulation
#'
#' Forwards model input to posterior simulation functions.
#'
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{create_df_model}}
#' in combination with \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of the class of the output of the applied posterior
#' simulation function. In case the package's own function is used, this will
#' result in an object of class \code{"dfm"}.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_df_model(x = bem_dfmdata, p = 1, n = 1,
#'                          iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#'
#' @export
draw_posterior.dfmodel <- function(object, ...){

  return(try(dfmpost(object)))
}
