#' Add Log-Likelihood
#' 
#' Generic function that calculates the posterior log-likelihoods of a model.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_posterior_loglik.bvarmodel}},
#' \code{\link{add_posterior_loglik.bvecmodel}},
#' \code{\link{add_posterior_loglik.expandingwindow}},
#' \code{\link{add_posterior_loglik.externalforecast}},
#' \code{\link{add_posterior_loglik.modellist}}.
#'
#' @export
add_posterior_loglik <- function (object, ...) {
  UseMethod("add_posterior_loglik")
}
