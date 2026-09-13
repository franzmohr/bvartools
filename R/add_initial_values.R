#' Add Initial Values of an MCMC Chain
#' 
#' A generic function used to generate initial values of an MCMC chain. The
#' function invokes particular methods which depend on the class of the first argument.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#' 
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_initial_values.bvarmodel}},
#' \code{\link{add_initial_values.bvecmodel}},
#' \code{\link{add_initial_values.expandingwindow}},
#' \code{\link{add_initial_values.externalforecast}},
#' \code{\link{add_initial_values.modellist}}.
#'
#' @export
add_initial_values <- function (object, ...) {
 UseMethod("add_initial_values")
}
