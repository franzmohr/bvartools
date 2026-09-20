#' @include multipliers.R
NULL

#' Dynamic Multipliers of a VEC Model with Exogenous Variables
#'
#' Computes the response of an endogenous variable to a change in a weakly
#' exogenous variable of an object of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel' with at least one exogenous
#' variable.
#' @param ... arguments of \code{\link{multipliers.bvarmodel}}, which does the
#' work.
#'
#' @details
#' A multiplier is a statement about the levels of the endogenous variables, so
#' the model is put into its levels form with \code{\link{vec_to_var}} and the
#' multipliers of that form are returned. The transformation is exact, and it is
#' where the error correction term does its work: the long-run relations enter
#' the levels coefficients, so a change in an exogenous variable that belongs to
#' a cointegrating relation moves the endogenous variables permanently, while at
#' rank zero it moves them only through the short-run terms.
#'
#' The responses are therefore in the units of the levels of the endogenous
#' variables, whatever the error correction form has on its left-hand side.
#'
#' @return A time-series object of class 'bvarirf', which is what
#' \code{\link{irf}} returns, so that the same \code{plot} method applies.
#'
#' @examples
#'
#' data("e6")
#' set.seed(1)
#' exogen <- ts(cbind(g = as.numeric(e6[, "R"]) * 0.3 + rnorm(nrow(e6), sd = 0.1)),
#'              start = start(e6), frequency = frequency(e6))
#'
#' model <- create_bvecmodel(data = e6, exogen = exogen, p = 2, s = 1, r = 1,
#'                           const = "unrestricted",
#'                           iterations = 100, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = 3, scale = 0.0001))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' dm <- multipliers(model, impulse = "g", response = "R", n_ahead = 8)
#'
#' @family post-estimation analysis
#' @export
#' @method multipliers bvecmodel
multipliers.bvecmodel <- function(x, ...) {

  if (is.null(x[["model"]][["m"]]) || x[["model"]][["m"]] == 0) {
    stop("Dynamic multipliers need a model with weakly exogenous variables, ",
         "and this model has none.")
  }
  if (is.null(x[["posterior"]][["a"]][["coeffs"]])) {
    stop("The model holds no posterior draws of its coefficients.")
  }

  multipliers(vec_to_var(x), ...)
}
