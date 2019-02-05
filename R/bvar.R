#' Bayesian Vector Autoregression Objects
#' 
#' `bvar` is used to create objects of class "bvar".
#' 
#' @param data the original time-series object of endogenous variables.
#' @param exogen the original time-series object of unmodelled variables.
#' @param y a \eqn{K \times T} matrix of endogenous variables,
#' usually, a result of a call to \code{\link{gen_var}}.
#' @param x a \eqn{(pK + (1+s)M + N) \times T} matrix of regressor variables, usually, a result of a
#' call to \code{\link{gen_var}}.
#' @param A0 a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters.
#' @param A a \eqn{pK^2 \times S} matrix of MCMC coefficient draws of lagged endogenous variables.
#' @param B a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables.
#' @param C an \eqn{KN \times S} matrix of MCMC coefficient draws of deterministic terms.
#' @param Sigma a \eqn{K^2 \times S} matrix of variance-covariance MCMC draws.
#' 
#' @details For the VARX model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_i y_{t-i} + \sum_{i = 0}^{s} B_i x_{t - i} + C D_t + u_t}
#' the function collects the S draws of a Gibbs sampler (after the burn-in phase) in a standardised object,
#' where \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of lagged endogenous variabels.
#' \eqn{x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{B_i} its corresponding coefficient matrix.
#' \eqn{D_t} is an N-dimensional vector of deterministic terms
#' and \eqn{C} its corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#' 
#' The draws of the different coefficient matrices provided in `A0`, `A`, `B`,
#' `C` and `Sigma` have to correspond to the same MCMC iteration.
#' 
#' @return An object of class "bvar" containing the following components, if specified:
#' \item{data}{the original time-series object of endogenous variables.}
#' \item{exogen}{the original time-series object of unmodelled variables.}
#' \item{y}{a \eqn{K \times T} matrix of endogenous variables.}
#' \item{x}{a \eqn{(pK + (1+s)M + N) \times T} matrix of regressor variables.}
#' \item{A0}{an \eqn{S \times K^2} `mcmc` object of coefficient draws of structural parameters.}
#' \item{A}{an \eqn{S \times pK^2} `mcmc` object of coefficient draws of lagged endogenous variables.}
#' \item{B}{an \eqn{S \times ((1 + s)MK)} `mcmc` object of coefficient draws of unmodelled, non-deterministic variables.}
#' \item{C}{an \eqn{S \times NK} `mcmc` object of coefficient draws of deterministic terms.}
#' \item{Sigma}{an \eqn{S \times K^2} `mcmc` object of variance-covariance draws.}
#' \item{specifications}{a list containing information on the model specification.}
#' 
#' @export
bvar <- function(data = NULL, exogen = NULL, y = NULL, x = NULL,
                 A0 = NULL, A = NULL, B = NULL,
                 C = NULL, Sigma = NULL) {

  result <- NULL
  if (is.null(y) | is.null(A)) {
    stop("At least the arguments 'y' and 'A' must be specified.")
  }
  if(!is.null(y)) {
    result$y <- y
    k <- NROW(y)
  }
  if(!is.null(A)) {
    result$A <- coda::mcmc(t(A))
    n_a <- ncol(result$A)
    if (n_a %% k == 0) {
      p <- n_a / k^2 
    } else {
      stop("Row number of argument 'A' is not a multiple of the number of endogenous variables.")
    }
  }
  if(!is.null(data)) {
    result$data <- data
  }
  if(!is.null(exogen)) {
    result$exogen <- exogen
  }
  if(!is.null(x)) {
    result$x <- x
  }
  if(!is.null(A0)) {
    result$A0 <- coda::mcmc(t(A0))
  }
  if(!is.null(B)) {
    result$B <- coda::mcmc(t(B))
  }
  if(!is.null(C)) {
    result$C <- coda::mcmc(t(C))
  }
  if(!is.null(Sigma)) {
    result$Sigma <- coda::mcmc(t(Sigma))
  }
  result$specifications <- list("dims" = c("K" = k),
                                "lags" = c("p" = p))
  
  class(result) <- append("bvar", class(result))
  return(result)
}