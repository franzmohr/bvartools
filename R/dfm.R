#' Bayesian Dynamic Factor Model Objects
#' 
#' \code{dfm} is used to create objects of class \code{"dfm"}.
#' 
#' @param x the standardised time-series object of observable variables.
#' @param lambda an \eqn{MN \times S} matrix of MCMC coefficient draws of factor loadings of the measurement equation.
#' @param fac an \eqn{NT \times S} matrix of MCMC draws of the factors in the transition equation, where the first N
#' rows correspond to the N factors in period 1 and the next N rows to the factors in period 2 etc.
#' @param sigma_u an \eqn{M \times S} matrix of MCMC draws for the error variances of the measurement equation.
#' @param a a \eqn{pN^2 \times S} matrix of MCMC coefficient draws of the transition equation.
#' @param sigma_v an \eqn{N \times S} matrix of MCMC draws for the error variances of the transition equation.
#' 
#' @details The function produces a standardised object from S draws of a Gibbs sampler (after the burn-in phase)
#' for the dynamic factor model (DFM) with measurement equation
#' \deqn{x_t = \lambda f_t + u_t,}
#' where
#' \eqn{x_t} is an \eqn{M \times 1} vector of observed variables,
#' \eqn{f_t} is an \eqn{N \times 1} vector of unobserved factors and
#' \eqn{\lambda} is the corresponding \eqn{M \times N} matrix of factor loadings.
#' \eqn{u_t} is an \eqn{M \times 1} error term.
#' 
#' The transition equation is
#' \deqn{f_t = \sum_{i=1}^{p} A_i f_{t - i} + v_t,}
#' where
#' \eqn{A_i} is an \eqn{N \times N} coefficient matrix and
#' \eqn{v_t} is an \eqn{N \times 1} error term.
#' 
#' @return An object of class \code{"dfm"} containing the following components, if specified:
#' \item{x}{the standardised time-series object of observable variables.}
#' \item{lambda}{an \eqn{S \times MN} "mcmc" object of draws of factor loadings of the measurement equation.}
#' \item{factor}{an \eqn{S \times NT} "mcmc" object of draws of factors.}
#' \item{sigma_u}{an \eqn{S \times M} "mcmc" object of variance draws of the measurement equation.}
#' \item{a}{an \eqn{S \times pN^2} "mcmc" object of coefficient draws of the transition equation.}
#' \item{sigma_v}{an \eqn{S \times N} "mcmc" object of variance draws of the transition equation.}
#' \item{specifications}{a list containing information on the model specification.}
#' 
#' @examples
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     lambda = list(v_i = .01),
#'                     sigma_u = list(shape = 5, rate = 4),
#'                     a = list(v_i = .01),
#'                     sigma_v = list(shape = 5, rate = 4))
#' 
#' # Obtain posterior draws
#' object <- dfmpost(model)
#' 
#' @export
dfm <- function(x, lambda = NULL, fac, sigma_u = NULL,
                a = NULL, sigma_v = NULL) {

  result <- NULL
  
  result$x <- x
  tt <- nrow(x)
  n <- nrow(fac) / tt
  
  if(!is.null(lambda)) {
    result$lambda <- coda::mcmc(t(lambda)) 
    n_lambda <- ncol(result$lambda)
  }
  
  result$factor <- coda::mcmc(t(fac))
  
  if(!is.null(sigma_u)) {
    result$sigma_u <- coda::mcmc(t(sigma_u)) 
  }
  
  p <- 0
  if(!is.null(a)) {
    result$a <- coda::mcmc(t(a))
    p <- NROW(a) / n^2
  }

  if(!is.null(sigma_v)) {
    result$sigma_v <- coda::mcmc(t(sigma_v)) 
  }
  
  result$specifications <- list("dims" = c("M" = NCOL(x), "N" = n),
                                "lags" = c("p" = p))
  
  class(result) <- "dfm"
  return(result)
}
