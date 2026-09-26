#' Posterior Simulation of Model Coefficients
#' 
#' Generic function used for posterior simulation of model coefficients.
#' 
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#'
#' @details Before the method is called, the size the object will have once its
#' draws are complete is compared with \code{options(bvartools.size_warning)},
#' and a warning is given if it is larger. See \code{\link{expected_size}}.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{add_posterior_coefficients.bvarmodel}},
#' \code{\link{add_posterior_coefficients.bvecmodel}},
#' \code{\link{add_posterior_coefficients.expandingwindow}},
#' \code{\link{add_posterior_coefficients.externalforecast}},
#' \code{\link{add_posterior_coefficients.modellist}}.
#'
#' @export
add_posterior_coefficients <- function(object, ...){
  # Once per call from the outside, for the whole of a list of models rather
  # than for each member the collection methods pass back through here.
  if (!isTRUE(.size_check[["posterior"]]) && .check_posterior_size(object, list(...)[["chains"]])) {
    .size_check[["posterior"]] <- TRUE
    on.exit(.size_check[["posterior"]] <- FALSE, add = TRUE)
  }
  UseMethod("add_posterior_coefficients")
}
