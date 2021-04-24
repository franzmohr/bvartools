#' Add Priors to Dynamic Factor Model
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{gen_dfm}}.
#'
#' @param object a list, usually, the output of a call to \code{\link{gen_dfm}}.
#' @param lambda a named list of prior specifications for the factor loadings in the measurement equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param sigma_u a named list of prior specifications for the error variance-covariance matrix. See 'Details'.
#' @param a a named list of prior specifications for the coefficients of the transition equation.
#' For the default specification the diagonal elements of the inverse prior variance-covariance matrix are set to 0.01.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' @param sigma_v a named list of prior specifications for the error variance-covariance matrix. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' Argument \code{lambda} can only contain the element \code{v_i}, which is a numeric specifying the prior
#' precision of the loading factors of the measurement equation. Default is 0.01.
#' 
#' The function assumes an inverse gamma prior for the errors of the measurement equation.
#' Argument \code{sigma_u} can contain the following elements:
#' \describe{
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error terms
#'   of the measurement equation. Default is 5.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error terms of the measurement
#'   equation. Default is 4.}
#' }
#' 
#' Argument \code{a} can only contain the element \code{v_i}, which is a numeric specifying the prior
#' precision of the coefficients of the transition equation. Default is 0.01.
#' 
#' The function assumes an inverse gamma prior for the errors of the transition equation.
#' Argument \code{sigma_v} can contain the following elements:
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
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples 
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 5000, burnin = 1000)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(v_i = .01),
#'                     sigma_u = list(shape = 5, rate = 4),
#'                     a = list(v_i = .01),
#'                     sigma_v = list(shape = 5, rate = 4))
#' 
#' @export
add_priors.dfmodel <- function(object,
                               lambda = list(v_i = 0.01),
                               sigma_u = list(shape = 5, rate = 4),
                               a = list(v_i = 0.01),
                               sigma_v = list(shape = 5, rate = 4),
                               ...){
  
  # rm(list = ls()[-which(ls() == "object")]); lambda = list(v_i = 0.01); sigma_v = list(shape = 1, rate = 1); phi = list(v_i = 0.01); sigma_u = list(shape = 1, rate = 1)
  
  only_one_model <- FALSE
  # If only one model is provided, make it compatible with the rest
  if ("data" %in% names(object)) {
    object <- list(object)
    only_one_model <- TRUE
  }
  
  # Checks - Coefficient priors ----
  if (!is.null(lambda)) {
    if (!is.null(lambda$v_i)) {
      if (lambda$v_i < 0) {
        stop("Argument 'lambda$v_i' must be at least 0.")
      }  
    } else {
      stop("Argument 'lambda$v_i' is missing.")
    }
  }
  
  if (!is.null(a)) {
    if (!is.null(a$v_i)) {
      if (a$v_i < 0) {
        stop("Argument 'a$v_i' must be at least 0.")
      }  
    } else {
      stop("Argument 'a$v_i' is missing.")
    }
  }
  
  # Checks - Error priors ----
  if (length(sigma_u) < 2) {
    stop("Argument 'sigma_u' must be at least of length 2.")
  } else {
    if (!"shape" %in% names(sigma_u)) {
      stop("Argument sigma_u$shape is missing.")
    }
    if (!"rate" %in% names(sigma_u)) {
      stop("Argument sigma_u$rate is missing.")
    }
    if (sigma_u$shape < 0) {
      stop("Argument 'sigma_u$shape' must be at least 0.")
    }
    if (sigma_u$rate <= 0) {
      stop("Argument 'sigma_u$rate' must be larger than 0.")
    } 
  }
  
  if (length(sigma_v) < 2) {
    stop("Argument 'sigma_v' must be at least of length 2.")
  } else {
    if (!"shape" %in% names(sigma_v)) {
      stop("Argument sigma_v$shape is missing.")
    }
    if (!"rate" %in% names(sigma_v)) {
      stop("Argument sigma_v$rate is missing.")
    }
    if (sigma_v$shape < 0) {
      stop("Argument 'sigma_v$shape' must be at least 0.")
    }
    if (sigma_v$rate <= 0) {
      stop("Argument 'sigma_v$rate' must be larger than 0.")
    } 
  }
  
  # Generate priors for each model specification ----
  for (i in 1:length(object)) {
    
    # Get model specs to obtain total number of coeffs
    m <- length(object[[i]]$model$factor$variables)
    n <- object[[i]]$model$factor$number
    p <- object[[i]]$model$factor$lags
    
    # Total number of freely estimated coefficients in lambda
    n_lambda <- (2 * m - n - 1) * n / 2
    
    # Total # of estimated coefficients in measurement equation
    n_a <- n * n * p
    
    # Priors for lambda ----
    object[[i]]$priors$lambda <- list(v_i = diag(lambda$v_i, n_lambda))
    
    # Priors for Phi ----
    if (n_a > 0) {
      object[[i]]$priors$a <- list(mu = matrix(0, n_a),
                                     v_i = diag(a$v_i, n_a))
    }
    
    # Error terms ----
    object[[i]]$priors$sigma_u$shape <- matrix(sigma_u$shape, m)
    object[[i]]$priors$sigma_u$rate <- matrix(sigma_u$rate, m)
    
    object[[i]]$priors$sigma_v$shape <- matrix(sigma_v$shape, n)
    object[[i]]$priors$sigma_v$rate <- matrix(sigma_v$rate, n)
  }
  
  if (only_one_model) {
    object <- object[[1]]
  }
  
  return(object)
}