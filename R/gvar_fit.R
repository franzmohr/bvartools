#' Fitter Function for GVAR Models
#' 
#' Estimates the country models of a GVAR model, combines them and solves the global model.
#' 
#' @param data a list of data and model specifications for country-specific models and information on the global
#' VAR model as produced by the \code{country_models} function.
#' @param iterations an integer of MCMC draws including burn-in (defaults to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler (defaults to 5000).
#' These draws do not enter the computation of posterior moments, forecasts etc.
#' @param thin an integer of thinning factor for MCMC output. Defaults to 10, which means that the
#' forecast sequences contain only every tenth draw of the original sequence. Set \code{thin = 1} to obtain the full MCMC sequence.
#' 
#' @details
#' The function applies the function \code{estimate_country_models} to each country model.
#' If \code{sfInit} was specified before, the estmation of each country model is distributed
#' across multiple cores by the load balancing function \code{sfClusterApplyLB}. Otherwise,
#' the estimations will be performed sequentially on a single core.
#' 
#' @return A list of data, model specifications, priors, thinned coefficient draws and information
#' criteria for each estimated country model, which are the output of the function \code{\link{estimate_country_model}}.
#' 
#' @export
gvar_fit <- function(data, iterations = 50000, burnin = 5000, thin = 10){
  
  if (!requireNamespace("snowfall")) {stop("The function 'gvar.fit' requires the packages 'snow' and 'snowfall'.")}
  
  #lag.max <- data[[1]]$specs$lags$maximum.lag
  #type <- data[[1]]$specs$type
  #case <- data[[1]]$specs$deterministic.terms
  #hrinkage <- data[[1]]$priors$Shrinkage[[1]]
  
  # Print estimation information
  cat(paste("Estimating a bunch of GVAR models.\n"))# with p = ", lag.max,
            #" and ", if (structural){"structural"}else{"reduced form"}, " ", type, " country models with\n",
            #if(tvp){"time varying"} else {"constant"}, " parameters, ",
            #if(sv){"stochasitc"} else {"constant"}, " volatility, deterministic case ", case, "\nand ",
            #if(shrinkage == "none"){"no"} else {paste(shrinkage, "as", sep = " ")}, " additional shrinkage method.\n", sep = ""))
  
  result <- snowfall::sfClusterApplyLB(data, estimate_country_model, iterations = iterations, burnin = burnin, thin = thin)
  names(result) <- names(data)
  
  return(result)
}