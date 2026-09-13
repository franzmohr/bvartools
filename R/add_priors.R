#' Add Priors to Bayesian Models
#'
#' A generic function used to generate prior vectors and matrices. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_priors.bvarmodel}},
#' \code{\link{add_priors.bvecmodel}}, \code{\link{add_priors.expandingwindow}},
#' \code{\link{add_priors.externalforecast}}, \code{\link{add_priors.modellist}}.
#'
#' @export
add_priors <- function (object, ...) {
 UseMethod("add_priors")
}
