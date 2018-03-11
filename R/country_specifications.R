#' Country Model Specifications
#' 
#' Produces a list of standard specifications for every country.
#' 
#' @param country.data a named list of "zoo" objects containing the country data.
#' @param global.data a "zoo" object with global data.
#' @param countries a character vector specifying the used countries in \code{country.data}.
#' @param p an integer or vector specifiying the lag order of domestic variables.
#' @param q an integer or vector specifiying the lag order of foreign variables.
#' @param s an integer or vector specifiying the lag order of global variables.
#' @param case a character specifying the deterministic terms used in the country models, see details below.
#' @param type a character specifying if the country models are estimated as vector autoregressive "VAR" (default) or error correction "VEC" models.
#' @param use a character vector specifying the considered variables. If \code{NULL}, all available variables will be used.
#' @param r an integer or vector for the cointegration rank, see details below.
#' @param structural a logical specifying if the country models should be estimated as structual models a la Primicieri (2005).
#' @param tvp a logical specifying if the country models should be estimated with time-varying parameters 
#' using the algorithm of Durbin & Koopman (2002).
#' @param sv a logical specifiying if the country models should be estimated with stochastic volatility
#' using the algorithm of Kaster & Frühwirth-Schnatter (2014).
#' 
#' @details
#' \code{p}, \code{q} and \code{s} refer to the lag order of the VAR formulation. If \code{type = "VEC"},
#' the function \code{\link{country.models}} will automatically reduce each lag order by 1. If a vector is provided,
#' \code{\link{country.models}} will produce a distinct model for all possible specification.
#' 
#' The argument \code{case} refers to Johansen (1995), where case I (default) refers to no deterministic terms and
#' case II and III specify a restricted and unrestricted constant term in the VEC formulation, respectively.
#' Both case II and case III enter a VAR model in the same way as an intercept. Case IV refers to a restricted trend and
#' an unrestricted constant in a VEC formulation. Case V adds an unrestricted constant and trend to the VEC model.
#' Again, cases IV and V are equivalent in a VAR model.
#' 
#' The argument \code{rank} is only used when \code{type = "VEC"}. Otherwise, the rank is set to zero.
#' If a vector is provided, \code{\link{country.models}} will produce a distinct model for all possible specifications.
#' 
#' @references
#' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for state space time series analysis. \emph{Biometrika}, 89(3), 603--615.
#'
#' Johansen, S. (1995) \emph{Likelihood-based inference in cointegrated vector autoregressive models}. Oxford University Press.
#'
#' Kastner, G. & Frühwirth-Schnatter, S. (2014). Ancillarity-sufficiency interweaving strategy (ASIS) for boosting MCMC estimation of stochastic volatility models. \emph{Computational Statistics & Data Analysis}, 76, 408–423, http://dx.doi.org/10.1016/j.csda.2013.01.002.
#' 
#' Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. \emph{The Review of Economic Studies}, 72(3), 821--852.
#' 
#' @return The function produces a list of preliminary model specifications, which consists of the following elements:
#' \item{domestic.variables}{a character vector of domestic variables in the country model.}
#' \item{foreign.variables}{a character vector of foreign variables in the country model.}
#' \item{global.variables}{a character vector of global variables in the country model.}
#' \item{deterministic.terms}{a character specifying the case referring to the used deterministic terms.}
#' \item{lags}{a list specifying the lag order of domestic, foreign and global variables and the order selection criterion.}
#' \item{rank}{a list specifying the cointegration rank of the country model.}
#' \item{type}{a character specifying the type of the model, either "VAR" or "VEC".}
#' 
#' @export
country_specifications <- function(country.data, global.data = NULL, countries = NULL,
                                   p = 2, q = 2, s = 2, case = "I", type = "VAR",
                                   use = NULL, r = NULL, structural = FALSE, tvp = FALSE, sv = FALSE){
  s.vars <- unique(unlist(lapply(country.data, names)))
  if (!is.null(use)){
    s.vars <- s.vars[which(is.element(s.vars, use))]
  }
  if (is.null(global.data)){
    g.vars <- NA
  } else {
    g.vars <- names(global.data)
    g.used <- which(is.element(g.vars, use))
    if (length(g.used) == 0){
      g.vars <- NA
      s <- NA
    } else {
      g.vars <- g.vars[g.used]
    }
  }
  if (type == "VAR"){
    r <- 0
  } else {
    if (is.null(r)){
      r <- NA 
    }
  }
  specs <- list()
  specs.names <- NULL
  for (i in 1:length(country.data)){
    d.vars <- names(country.data[[i]])
    if (!is.null(use)){
      d.vars <- d.vars[which(is.element(d.vars, use))]
    }
    specs <- c(specs, list(list(domestic.variables = d.vars,
                                foreign.variables = s.vars,
                                global.variables = g.vars,
                                deterministic.terms = case,
                                lags = list(lags = list("domestic" = p, "foreign" = q, "global" = s),
                                            fixed = list("domestic" = FALSE, "foreign" = FALSE, "global" = FALSE)),
                                rank = list(rank = r),
                                type = type,
                                structural = structural,
                                tvp = tvp,
                                sv = sv)))
    specs.names <- c(specs.names, names(country.data)[i])
  }
  names(specs) <- specs.names
  for (i in 1:length(specs)){
    specs[[i]][[1]] <- specs[[i]][[1]][is.element(specs[[i]][[1]], names(country.data[[i]]))]
  }
  if (!is.null(countries)){
    if (any(!is.element(countries, names(country.data)))) {
      stop(paste("There are no country data available for ", paste(countries[which(!is.element(countries, names(country.data)))], collapse = ", "), ".", sep = ""))
    }
    result <- NULL
    names.result <- NULL
    for (i in countries){
      result <- c(result, list(specs[[i]]))
      names.result <- c(names.result, i)
    }
    names(result) <- names.result
  } else {
    result <- specs
  }
  return(result)
}