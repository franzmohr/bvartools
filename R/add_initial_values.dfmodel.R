#' Add Initial Values to a Dynamic Factor Model
#'
#' Adds initial values to a dynamic factor model, which was produced by
#' function \code{\link{create_df_model}} in combination with \code{\link{add_priors}}.
#'
#' @param object a named list, usually, the output of a call to \code{\link{create_df_model}}.
#' @param method a character specifying the method of how initial values are generated.
#' Defaults to \code{"prior"}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' For argument \code{method} the following specifications are possible:
#' \describe{
#'   \item{\code{"prior"}}{Initial values are drawn from the prior. Not possible for uninformative priors.}
#' }
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_df_model(x = bem_dfmdata, p = 1:2, n = 1,
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
#' @export
add_initial_values.dfmodel <- function(object, method = "prior", ...){
  
  if (method == "prior") {
    
    # lambda
    object$initial$lambda <- chol(object$priors$lambda$vinv) %*% stats::rnorm(nrow(object$priors$lambda$vinv))
    
    # U
    sigma_shape <- object$priors$u$shape
    sigma_rate <- 1 / object$priors$u$rate
    object$initial$uinv <- diag(1, object$model$m)
    for (i in 1:object$model$m) {
      object$initial$uinv[i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
    }
    rm(list = c("sigma_shape", "sigma_rate"))
    
    # V
    sigma_shape <- object$priors$v$shape
    sigma_rate <- 1 / object$priors$v$rate
    object$initial$vinv <- diag(1, object$model$n)
    for (i in 1:object$model$n) {
      object$initial$vinv[i, i] <- 1 / stats::rgamma(1, shape = sigma_shape[i], rate = sigma_rate[i])
    }
    rm(list = c("sigma_shape", "sigma_rate"))
    
    if (object$model$p > 0) {
      # A
      object$initial$a <- object$priors$a$mu + chol(object$priors$a$vinv) %*% stats::rnorm(object$model$n^2 * object$model$p)
    }
  }
  
  return(object)
}
