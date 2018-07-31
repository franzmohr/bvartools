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
#' 
#' @return A list of data, model specifications, priors, thinned coefficient draws and the values of information criteria.
#' 
#' @export
estimate_country_model <- function(data, iterations = 50000, burnin = 5000, thin = 10) {
  if (type == "VAR"){
    stop("BVAR not implemented yet.")
    result <- try(bvarx(data, iterations = iterations, burnin = burnin, thin = thin))
  }
  if (type == "VEC") {
    if (is.na(data$specs$rank)) {
      result <- try(bvecx(data, iterations = iterations, burnin = burnin, thin = thin)) 
    } else {
      stop("Reduced rank BVEC not implemented yet.")
      #result <- try(bvecx_rr(data, iterations = iterations, burnin = burnin, thin = thin)) 
    }
  }
  if (inherits(result, "try-error")) {
    result <- c(data, list(coefs = NULL, criteria = NA))
    }
  return(result)
}