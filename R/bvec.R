#' Bayesian Vector Error Correction Objects
#' 
#' The function `bvec` is used to create objects of class "bvec".
#'
#' @param y a \eqn{k x T} matrix of differenced dependent variables.
#' @param ect a matrix of cointegrating variables.
#' @param x a matrix of differenced regressors.
#' @param A0 a \eqn{k^2 x S} matrix of structural coefficients, where \eqn{S} is the
#' total amount of MCMC draws.
#' @param alpha a matrix of MCMC draws.
#' @param beta a matrix of MCMC draws.
#' @param Pi a matrix of MCMC draws.
#' @param Gamma a matrix of MCMC draws.
#' @param Ypsilon a matrix of MCMC draws.
#' @param D a matrix of MCMC draws.
#' @param Sigma a matrix of MCMC draws.
#' @param LL a matrix of MCMC draws.
#' 
#' @details The function collects the output of a Gibbs sampler in a standardised object,
#' which can be used as input for further analyis such as forecasting and impulse
#' response analysis.
#' 
#' The notation corresponds to the following model
#' \deqn{y_t = \Pi ECT_{t} +  \sum_{i=1}^p Gamma_i y_{t-i} + \sum_{i=0}^s Ypsilon_i x_{t-i} + A_0 u_t,}
#' with \eqn{\Pi = \alpha \beta'} and \eqn{u_t \sim \Sigma}.
#' 
#' @return A list containing the following elements:
#' \item{y}{a \eqn{n x T} matrix.}
#' \item{ect}{error correction term.}
#' \item{x}{a \eqn{pn + qm x T} matrix.}
#' \item{A0}{a \eqn{D x pn^2} `mcmc` object of structural parameters draws.}
#' \item{alpha}{a `mcmc` object of loading matrix.}
#' \item{beta}{a `mcmc` object of cointegration matrix.}
#' \item{Pi}{a `mcmc` object of cointegration matrix.}
#' \item{Gamma}{a \eqn{D x pn^2} `mcmc` object of lagged AR coefficient draws.}
#' \item{Ypsilon}{a \eqn{D x pn^2} `mcmc` object of contemporanesous and lagged exogenous coefficient draws.}
#' \item{D}{a \eqn{D x dn} `mcmc` object of deterministic terms.}
#' \item{Sigma}{a \eqn{D x n^2} `mcmc` object of variance-covariance draws.}
#' \item{LL}{a \eqn{D x t} `mcmc` object of log-likelihood draws.}
#' 
#' @export
bvec <- function(y = NULL, ect = NULL, x = NULL, A0 = NULL,
                 alpha = NULL, beta = NULL, Pi = NULL,
                 Gamma = NULL, Ypsilon = NULL, D = NULL, Sigma = NULL, LL = NULL) {

  result <- NULL
  
  if(!is.null(y)) {
    result$y <- y
  }
  if(!is.null(ect)) {
    result$ect <- ect
  }
  if(!is.null(x)) {
    result$x <- x
  }
  if(!is.null(A0)) {
    result$A0 <- coda::mcmc(t(A0))
  }
  if(!is.null(alpha)) {
    result$alpha <- coda::mcmc(t(alpha))
  }
  if(!is.null(beta)) {
    result$beta <- coda::mcmc(t(beta))
  }
  if(!is.null(Pi)) {
    result$Pi <- coda::mcmc(t(Pi))
  }
  if(!is.null(Gamma)) {
    result$Gamma <- coda::mcmc(t(Gamma))
  }
  if(!is.null(Ypsilon)) {
    result$Ypsilon <- coda::mcmc(t(Ypsilon))
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
  
  class(result) <- append("bvec", class(result))
  return(result)
}