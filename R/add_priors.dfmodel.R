#' Add Priors to Bayesian Models
#' 
#' Adds prior specifications to a dynamic factor model, which was produced by
#' function \code{\link{create_df_model}}.
#'
#' @param object a list of class 'dfmodel'.
#' @param lambda a named list of prior specifications for the factor loadings in the measurement equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param u a named list of prior specifications for the error variance-covariance matrix. See 'Details'.
#' @param a a named list of prior specifications for the coefficients of the transition equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param v a named list of prior specifications for the error variance-covariance matrix. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Argument \code{lambda} can only contain the element \code{vinv}, which is a numeric specifying the prior
#' precision of the loading factors of the measurement equation. Default is 0.01.
#'
#' The function assumes an inverse gamma prior for the errors of the measurement equation.
#' Argument \code{u} can contain the following elements:
#' \describe{
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error terms
#'   of the measurement equation. Default is 5.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error terms of the measurement
#'   equation. Default is 4.}
#' }
#'
#' Argument \code{a} can only contain the element \code{vinv}, which is a numeric specifying the prior
#' precision of the coefficients of the transition equation. Default is 0.01.
#'
#' The function assumes an inverse gamma prior for the errors of the transition equation.
#' Argument \code{v} can contain the following elements:
#' \describe{
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error terms
#'   of the transition equation. Default is 5.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error terms of the transition
#'   equation. Default is 4.}
#' }
#'
#' @return A list of models.
#'
#' @references
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#'
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' @examples
#'
#' # Load data
#' data("bem_dfmdata")
#'
#' # Generate model data
#' model <- create_df_model(x = bem_dfmdata, p = 1:2, n = 1,
#'                          iterations = 5000, burnin = 1000)
#' # Number of iterations and burn-in should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(vinv = .01),
#'                     u = list(shape = 5, rate = 4),
#'                     a = list(vinv = .01),
#'                     v = list(shape = 5, rate = 4))
#'
#' @export
add_priors.dfmodel <- function(object,
                               lambda = list(vinv = 0.01),
                               u = list(shape = 5, rate = 4),
                               a = list(vinv = 0.01),
                               v = list(shape = 5, rate = 4),
                               ...){

  # Checks - Coefficient priors ----
  if (!is.null(lambda)) {
    if (!is.null(lambda$vinv)) {
      if (lambda$vinv < 0) {
        stop("Argument 'lambda$vinv' must be at least 0.")
      }
    } else {
      stop("Argument 'lambda$vinv' is missing.")
    }
  }

  if (!is.null(a)) {
    if (!is.null(a$vinv)) {
      if (a$vinv < 0) {
        stop("Argument 'a$vinv' must be at least 0.")
      }
    } else {
      stop("Argument 'a$vinv' is missing.")
    }
  }

  # Checks - Error priors ----
  if (length(u) < 2) {
    stop("Argument 'u' must be at least of length 2.")
  } else {
    if (!"shape" %in% names(u)) {
      stop("Argument u$shape is missing.")
    }
    if (!"rate" %in% names(u)) {
      stop("Argument u$rate is missing.")
    }
    if (u$shape < 0) {
      stop("Argument 'u$shape' must be at least 0.")
    }
    if (u$rate <= 0) {
      stop("Argument 'u$rate' must be larger than 0.")
    }
  }

  if (length(v) < 2) {
    stop("Argument 'v' must be at least of length 2.")
  } else {
    if (!"shape" %in% names(v)) {
      stop("Argument v$shape is missing.")
    }
    if (!"rate" %in% names(v)) {
      stop("Argument v$rate is missing.")
    }
    if (v$shape < 0) {
      stop("Argument 'v$shape' must be at least 0.")
    }
    if (v$rate <= 0) {
      stop("Argument 'v$rate' must be larger than 0.")
    }
  }

  # Get model specs to obtain total number of coeffs
  m <- length(object$model$variables)
  n <- object$model$n_factors
  p <- object$model$p

  # Total number of freely estimated coefficients in lambda
  n_lambda <- (2 * m - n - 1) * n / 2

  # Total # of estimated coefficients in measurement equation
  n_a <- n * n * p

  # Priors for lambda ----
  object$priors$lambda <- list(vinv = diag(lambda$vinv, n_lambda))

  # Priors for Phi ----
  if (n_a > 0) {
    object$priors$a <- list(mu = matrix(0, n_a),
                                 vinv = diag(a$vinv, n_a))
  }

  # Error terms ----
  object$priors$u$shape <- matrix(u$shape, m)
  object$priors$u$rate <- matrix(u$rate, m)

  object$priors$v$shape <- matrix(v$shape, n)
  object$priors$v$rate <- matrix(v$rate, n)

  return(object)
}
