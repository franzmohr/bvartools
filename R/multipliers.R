#' Dynamic Multipliers
#'
#' A generic function used to calculate the dynamic multipliers of a model with
#' weakly exogenous variables. The function invokes particular methods which
#' depend on the class of the first argument.
#'
#' @param x an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @details
#' A dynamic multiplier is the response of the endogenous variables of a model
#' to a change in one of its weakly exogenous variables. It is what an impulse
#' response is for a shock to an error term, for a variable the model does not
#' explain: the exogenous variable is moved by hand and the endogenous ones are
#' followed.
#'
#' The quantity is central to models whose foreign block is weakly exogenous,
#' such as the country models of a global VAR (Pesaran, Schuermann and Weiner,
#' 2004), where the responses to a change in the foreign variables are what
#' the country model has to say about the rest of the world without the rest
#' of the world being solved.
#'
#' @return The value returned by the method for the class of \code{x}, as
#' described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{multipliers.bvarmodel}},
#' \code{\link{multipliers.bvecmodel}}.
#'
#' @references
#'
#' Pesaran, M. H., Schuermann, T., Weiner, S. M. (2004). Modeling regional
#' interdependencies using a global error-correcting macroeconometric model.
#' \emph{Journal of Business & Economic Statistics, 22}(2), 129-162.
#'
#' @family post-estimation analysis
#' @export
multipliers <- function(x, ...) {
  UseMethod("multipliers")
}
