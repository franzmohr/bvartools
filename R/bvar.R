#' Bayesian Vector Autoregression Objects
#' 
#' The function `bvar` is used to create objects of class "bvar".
#'
#' @param y a \eqn{k x T} matrix, where \eqn{k} is the number of endogenous variables and
#' \eqn{T} is the total amount of observations.
#' @param x a \eqn{pk x T} matrix, where \eqn{p} is the lag order of the VAR model.
#' @param A0 a \eqn{k^2 x S} matrix of structural coefficients, where \eqn{S} is the
#' total amount of MCMC draws.
#' @param A a \eqn{pn^2 x S} matrix of lagged AR coefficients.
#' @param B a \eqn{(1 + q)mn x S} matrix of coefficients of contemporaneous and
#' lagged exogenous variables.
#' @param D a \eqn{dn x S} matrix, where \eqn{d} is the number of deterministic
#' variables.
# @param M a \eqn{n^2 x S} matrix.
#' @param Sigma a \eqn{n^2 x S} matrix of variance-covariance draws.
#' @param LL a \eqn{t x S} matrix of log-likelihood draws.
#' 
#' @details The function collects the output of a Gibbs sampler in a standardised object,
#' which can be used as input for further analyis such as forecasting and impulse
#' response analysis.
#' 
#' The notation corresponds to the following model
#' \deqn{y_t = \sum_{i=1}^p A_i y_{t-i} + \sum_{i=0}^q B_i x_{t-i} + A_0 u_t,}
#' with \eqn{u_t \sim \Sigma}.
#' 
#' @return A list containing the following elements:
#' \item{y}{a \eqn{n x T} matrix.}
#' \item{x}{a \eqn{pn + qm x T} matrix.}
#' \item{A0}{a \eqn{D x pn^2} `mcmc` object of structural parameters draws.}
#' \item{A}{a \eqn{D x pn^2} `mcmc` object of lagged AR coefficient draws.}
#' \item{B}{a \eqn{D x pn^2} `mcmc` object of contemporanesous and lagged exogenous coefficient draws.}
#' \item{D}{a \eqn{D x dn} `mcmc` object of deterministic terms.}
#' \item{Sigma}{a \eqn{D x n^2} `mcmc` object of variance-covariance draws.}
#' \item{LL}{a \eqn{D x t} `mcmc` object of log-likelihood draws.}
#' 
#' @export
bvar <- function(y = NULL, x = NULL, A0 = NULL, A = NULL, B = NULL,
                  D = NULL, Sigma = NULL, LL = NULL) {

  result <- NULL
  
  if(!is.null(y)) {
    result$y <- y
  }
  if(!is.null(x)) {
    result$x <- x
  }
  if(!is.null(A0)) {
    result$A0 <- coda::mcmc(t(A0))
  }
  if(!is.null(A)) {
    result$A <- coda::mcmc(t(A))
  }
  if(!is.null(B)) {
    result$B <- coda::mcmc(t(B))
  }
  if(!is.null(D)) {
    result$D <- coda::mcmc(t(D))
  }
  if(!is.null(Sigma)) {
    result$Sigma <- coda::mcmc(t(Sigma))
  }
  if(!is.null(LL)) {
    result$LL <- coda::mcmc(t(LL))
  }
  
  class(result) <- append("bvar", class(result))
  return(result)
}