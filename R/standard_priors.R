#' Standard Priors
#' 
#' Produces a list of standard priors used in the estimation of the country models.
#' 
#' @return A list containing a set of standard priors.
#' \item{A}{a list with priors for constant and time varying coefficients of non-deterministic variables.}
#' \item{Deterministic}{a list with priors for constant and time-varying coefficients of deterministic variables.}
#' \item{A0}{a list with priors for constant and time-varying coefficients of structural variables.}
#' \item{Pi}{a list with priors for constant and time-varying coefficients of variables in the cointegration term.}
#' \item{Sigma}{a list with priors for the estimation of constant covariance matrices and stochastic volatility.}
#' \item{Shrinkage}{a list with type and priors for additional shinkage methods.}
#' 
#' @export
standard_priors <- function(){
  result <- list("A" = list("constant" = c(0, 1 / 2), 
                            "tvp" = c(3, .0001)),
                 "Deterministic" = list("constant" = c(0, 1 / 9),
                                        "tvp" = c(3, .0001)),
                 "A0" = list("constant" = c(0, 1 / 2),
                             "tvp" = c(3, .0001)),
                 "Pi" = list("non_rr" = list("constant" = c(0, 1 / 2),
                                             "tvp" = c(3, .0001)),
                             "reduced_rank" = list("constant" = list("alpha" = 0, "v_i" = 0, "P_i" = 1, "G_i" = "Omega_i"),
                                                   "tvp" = list("alpha" = c(3, .0001),
                                                                "rho" = c(.999, 1)))),
                 "Omega" = list("constant" = c(1, 0.000001)),
                 "Sigma" = list("sv" = list(para = list(mu = -9,
                                                        phi = .95,
                                                        sigma = .1),
                                            latent = -9,
                                            priors = list(priormu = c(-9, 10), 
                                                          priorphi = c(25, 1.5),
                                                          priorsigma = .2))),
                 "Shrinkage" = list("type" = "none"))
  return(result)
}