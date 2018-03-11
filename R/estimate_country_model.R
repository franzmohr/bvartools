#' Estimate a Country Model
#' 
#' Assigns the framework, i.e., data, specifications and priors, of a country model to a specified estimator function. The
#' function ensures that the estimation process is not stopped when an estimator function produces an error.
#' 
#' @param data a list containing the data, weights, specifications and priors of a country model, i.e. a component
#' of the list produced by \code{\link{country_models}}.
#' @param iterations an integer of MCMC draws including burn-in (defaults to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler (defaults to 5000). These draws do
#' not enter the computation of posterior moments, forecasts etc.
#' @param thin an integer of thinning factor for MCMC output. Defaults to 10, which means that the forecast
#' sequences contain only every tenth draw of the original sequence. Set \code{thin = 1} to obtain the full MCMC sequence.
#' @param structural a logical specifying if the country models should be estimated as structual models a la Primicieri (2005).
#' @param tvp a logical specifying if the country models should be estimated with time-varying parameters
#' using the algorithm of Durbin & Koopman (2002).
#' @param sv a logical specifying if the country models should be estimated with stochastic volatility
#' using the algorithm of Kaster & Frühwirth-Schnatter (2014).
#' 
#' @return A list of data, model specifications, priors, thinned coefficient draws and the values of information criteria.
#' 
#' @export
estimate_country_model <- function(data, iterations = 50000, burnin = 5000, thin = 10) {
  type <- data$specs$type
  structural <- data$specs$structural
  tvp <- data$specs$tvp
  sv <- data$specs$sv
  if (type == "VAR"){
    if (structural){
      if (tvp){
        #result <- snowfall::sfClusterApplyLB(data, tvp.bsvarx, iterations = iterations, burnin = burnin, thin = thin, sv = sv)
      } else {
        #result <- snowfall::sfClusterApplyLB(data, bsvarx, iterations = iterations, burnin = burnin, thin = thin, sv = sv)
      }
    } else {
      if (tvp){
        #result <- snowfall::sfClusterApplyLB(data, tvp.bvarx, iterations = iterations, burnin = burnin, thin = thin)
      } else {
        result <- try(bvarx(data, iterations = iterations, burnin = burnin, thin = thin), silent = TRUE)
      }
    }
  }
  if (type == "VEC") {
    if (structural){
      if (tvp){
        #result <- snowfall::sfClusterApplyLB(data, tvp.bsvecx, iterations = iterations, burnin = burnin, thin = thin, sv = sv)
      } else {
        #result <- snowfall::sfClusterApplyLB(data, bsvecx, iterations = iterations, burnin = burnin, thin = thin, sv = sv)
      }
    } else {
      if (tvp){
        #result <- snowfall::sfClusterApplyLB(data, tvp.bvecx, iterations = iterations, burnin = burnin, thin = thin)
      } else {
        result <- bvecx(data, iterations = iterations, burnin = burnin, thin = thin)
      }
    }
  }
  
  if (inherits(result, "try-error")) {result <- c(data, list(coef = NULL, criteria = NA))}
  return(result)
}