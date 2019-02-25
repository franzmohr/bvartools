#' Bayesian Vector Error Correction Objects
#' 
#' `bvec` is used to create objects of class "bvec".
#' 
#' @param data the original time-series object of endogenous variables.
#' @param exogen the original time-series object of unmodelled variables.
#' @param y a \eqn{K \times T} matrix of differenced endogenous variables,
#' usually, a result of a call to \code{\link{gen_vec}}.
#' @param w a \eqn{(K + M + N^{R}) \times T} matrix of variables in the
#' cointegration term, usually, a result of a call to \code{\link{gen_vec}}.
#' @param x a \eqn{(K(p - 1) + Ms + N^{UR}) \times T} matrix of differenced regressors
#' of \eqn{y} and \eqn{x}, and unrestricted deterministic terms, usually,
#' a result of a call to \code{\link{gen_vec}}.
#' @param A0 a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters.
#' @param alpha a \eqn{Kr \times S} matrix of MCMC coefficient draws of the loading matrix \eqn{\alpha}.
#' @param beta a \eqn{((K + M + N^{R})r) \times S} matrix of MCMC coefficient draws of cointegration matrix \eqn{\beta}.
#' @param Pi a \eqn{K^2 \times S} matrix of MCMC coefficient draws of endogenous varaibles in the cointegration matrix.
#' @param Pi_x a \eqn{KM \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables in the cointegration matrix.
#' @param Pi_d a \eqn{KN^{R} \times S} matrix of MCMC coefficient draws of restricted deterministic terms.
#' @param Gamma a \eqn{(p-1)K^2 \times S} matrix of MCMC coefficient draws of differenced lagged endogenous variables.
#' @param Upsilon an \eqn{sMK \times S} matrix of MCMC coefficient draws of differenced unmodelled variables.
#' @param C an \eqn{KN^{UR} \times S} matrix of MCMC coefficient draws of unrestricted deterministic terms.
#' @param Sigma a \eqn{K^2 \times S} matrix of variance-covariance MCMC draws.
#' 
#' @details For the VECX model
#' \deqn{A_0 \Delta y_t = \Pi^{+} \begin{pmatrix} y_{t-1} \\ x_{t-1} \\ d^{R}_{t-1} \end{pmatrix} +
#' \sum_{i = 1}^{p-1} \Gamma_i \Delta y_{t-i} +
#' \sum_{i = 0}^{s-1} \Upsilon_i \Delta x_{t-i} +
#' C^{UR} d^{UR}_t + u_t}
#' the function collects the S draws of a Gibbs sampler (after the burn-in phase) in a standardised object,
#' where \eqn{\Delta y_t} is a K-dimensional vector of differenced endogenous variables
#' and \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{\Pi^{+} = \left[ \Pi, \Pi^{x}, \Pi^{d} \right]} is
#' the coefficient matrix of the error correction term, where
#' \eqn{y_{t-1}}, \eqn{x_{t-1}} and \eqn{d^{R}_{t-1}} are the first lags of endogenous,
#' exogenous variables in levels and restricted deterministic terms, respectively.
#' \eqn{\Pi}, \eqn{\Pi^{x}}, and \eqn{\Pi^{d}} are the corresponding coefficient matrices, respectively.
#' \eqn{\Gamma_i} is a coefficient matrix of lagged differenced endogenous variabels.
#' \eqn{\Delta x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{\Upsilon_i} its corresponding coefficient matrix. \eqn{d_t} is an
#' \eqn{N^{UR}}-dimensional vector of unrestricted deterministics and \eqn{C^{UR}}
#' the corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#' 
#' The draws of the different coefficient matrices provided in \code{alpha}, \code{beta},
#' \code{Pi}, \code{Pi_x}, \code{Pi_d}, \code{A0}, \code{Gamma}, \code{Ypsilon},
#' \code{C} and \code{Sigma} have to correspond to the same MCMC iteration.
#' 
#' @return An object of class "gvec" containing the following components, if specified:
#' \item{data}{the original time-series object of endogenous variables.}
#' \item{exogen}{the original time-series object of unmodelled variables.}
#' \item{y}{a \eqn{K \times T} matrix of differenced endogenous variables.}
#' \item{w}{a \eqn{(K + M + N^{R}) \times T} matrix of variables in the cointegration term.}
#' \item{x}{a \eqn{((p - 1)K + sM + N^{UR}) \times T} matrix of differenced regressor variables and
#' unrestricted deterministic terms.}
#' \item{A0}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of structural parameters.}
#' \item{alpha}{an \eqn{S \times Kr} "mcmc" object of coefficient draws of loading parameters.}
#' \item{beta}{an \eqn{S \times ((K + M + N^{R})r)} "mcmc" object of coefficient draws of cointegration parameters.}
#' \item{Pi}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of endogenous variables in the cointegration matrix.}
#' \item{Pi_x}{an \eqn{S \times KM} "mcmc" object of coefficient draws of unmodelled, non-deterministic variables in the cointegration matrix.}
#' \item{Pi_d}{an \eqn{S \times KN^{R}} "mcmc" object of coefficient draws of unrestricted deterministic variables in the cointegration matrix.}
#' \item{Gamma}{an \eqn{S \times (p-1)K^2} "mcmc" object of coefficient draws of differenced lagged endogenous variables.}
#' \item{Upsilon}{an \eqn{S \times sMK} "mcmc" object of coefficient draws of differenced unmodelled variables.}
#' \item{C}{an \eqn{S \times KN^{UR}} "mcmc" object of coefficient draws of deterministic terms.}
#' \item{Sigma}{an \eqn{S \times K^2} "mcmc" object of variance-covariance draws.}
#' \item{specifications}{a list containing information on the model specification.}
#' 
#' @export
bvec <- function(data = NULL, exogen = NULL, y = NULL, w = NULL, x = NULL,
                 alpha = NULL, beta = NULL, Pi = NULL, Pi_x = NULL, Pi_d = NULL,
                 A0 = NULL, Gamma = NULL, Upsilon = NULL, C = NULL, Sigma = NULL) {

  result <- NULL
  if (is.null(y) | is.null(Pi)) {
    stop("At least the arguments 'y' and 'Pi' must be specified.")
  }
  if(!is.null(y)) {
    result$y <- y
    k <- NROW(y)
  }
  if(!is.null(Pi)) {
    result$Pi <- coda::mcmc(t(Pi))
    if (NROW(Pi) != k^2) {
      stop("Row number of argument 'Pi' is not equal to the number of endogenous variables.")
    }
  }
  if(!is.null(w)) {
    result$w <- w
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
  if(!is.null(Pi_x)) {
    result$Pi_x <- coda::mcmc(t(Pi_x))
  }
  if(!is.null(Pi_d)) {
    result$Pi_d <- coda::mcmc(t(Pi_d))
  }
  if(!is.null(Gamma)) {
    result$Gamma <- coda::mcmc(t(Gamma))
    n_gamma <- ncol(result$Gamma)
    if (n_gamma %% k == 0) {
      p <- n_gamma / k^2
    } else {
      stop("Row number of argument 'Gamma' is not a multiple of the number of endogenous variables.")
    }
  } else {
    p <- 0
  }
  if(!is.null(Upsilon)) {
    result$Ypsilon <- coda::mcmc(t(Upsilon))
  }
  if(!is.null(C)) {
    result$C <- coda::mcmc(t(C))
  }
  if(!is.null(Sigma)) {
    result$Sigma <- coda::mcmc(t(Sigma))
  }
  result$specifications <- list("dims" = c("K" = k),
                                "lags" = c("p" = p + 1))
  
  class(result) <- append("bvec", class(result))
  return(result)
}