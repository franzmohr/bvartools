#' Persistence Profiles
#'
#' A generic function used to calculate the persistence profiles of the
#' cointegrating relations of a model.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{persistence_profiles.bvecmodel}}.
#'
#' @references
#'
#' Pesaran, M. H., & Shin, Y. (1996). Cointegration and speed of convergence to
#' equilibrium. \emph{Journal of Econometrics, 71}(1-2), 117--143.
#' \doi{10.1016/0304-4076(94)01697-6}
#'
#' @export
persistence_profiles <- function(object, ...) {
  UseMethod("persistence_profiles")
}
