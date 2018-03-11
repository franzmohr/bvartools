#' Add Bayesian Variable Selection Priors to Country Models
#' 
#' Produces a 
#' 
#' @param data a list as produced by \link{\code{country_models}}.
#' @param method
add_BVS_prior <- function(data){
  x <- data$Country.data$USA
  
  if (x$specs$type == "VAR") {
    
  }
  
  return(data)
}