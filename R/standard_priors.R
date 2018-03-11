#' Standard Priors
#' 
#' Produces a list of standard priors used in the estimation of the country models.
#' 
#' @param shrinkage a character describing the shrinkage method. Either \code{NULL} (default) for none, "SSVS"
#' for stochastic search variable selection as in George et. al (2008) or "BVS" for Bayesian variable selection
#' as in Korobilis (2013).
#' 
#' @return A list containing a set of standard priors.
#' \item{A}{a list with priors for constant and time varying coefficients of non-deterministic variables.}
#' \item{Deterministic}{a list with priors for constant and time-varying coefficients of deterministic variables.}
#' \item{A0}{a list with priors for constant and time-varying coefficients of structural variables.}
#' \item{Pi}{a list with priors for constant and time-varying coefficients of variables in the cointegration term.}
#' \item{Sigma}{a list with priors for the estimation of constant covariance matrices and stochastic volatility.}
#' \item{Shrinkage}{a list with type and priors for additional shinkage methods.}
#' 
#' @references 
#' George, E. I.; Sun, D. & Ni, S. (2008). Bayesian stochastic search for VAR model restrictions. \emph{Journal of Econometrics}, 142(1), 553--580.
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection. \emph{Journal of Applied Econometrics}, 28(2), 204--230.
#' 
#' @export
standard_priors <- function(shrinkage = NULL){
  result <- list("A" = list("constant" = c(0, 1/9), 
                            "tvp" = c(2, .0001)),
                 "Deterministic" = list("constant" = c(0, 1/9),
                                        "tvp" = c(2, .0001)),
                 "A0" = list("constant" = c(0, 1/9),
                             "tvp" = c(2, .0001)),
                 "Pi" = list("alpha" = list("constant" = c(0, 0),
                                            "tvp" = c(2, .0001)),
                             "beta" = list("constant" = c(0, 0),
                                           "tvp" = c(2, .0001)),
                             "rho" = .9999),
                 "Sigma" = list("constant" = c(1, 0.0001),
                                "sv" = list(para = list(mu = -9, phi = .9, sigma = .1),
                                            latent = 0,
                                            priors = list(priormu = c(-9, 10), 
                                                          priorphi = c(25, 1.5),
                                                          priorsigma = .2))),
                 "Shrinkage" = list("type" = "none"))
  
  if (!is.null(shrinkage)){
    if (shrinkage == "BVS") {
      result$Shrinkage <- list("type" = "BVS", "spec" = .5, "exclude.deterministic" = TRUE)
    }
    if (shrinkage == "SSVS") {
      result$Shrinkage <- list("type" = "SSVS", "spec" = .5, "exclude.deterministic" = TRUE)
    }
  }
  
  return(result)
}