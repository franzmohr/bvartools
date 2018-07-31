#' Specifications of Country Models in a GVAR Model
#' 
#' Extracts the name, lags of domestic, foreign and global variables and the rank of the cointegration
#' matrix of all country models, which were used for the construction of the global model.
#' 
#' @param data a list containing the individual country model results and the solved
#' GVAR model as produced by \code{\link{solve_gvar}}.
#' 
#' @return A data frame containing the specification for each country model,
#' which was used to construct the global model.
#' 
#' @export
get_country_specs <- function(data){
  result <- c()
  for (i in which(data$GVAR$specs$criteria[, "Maximum"])) {
    result <- rbind(result,
                    data.frame("Country" = names(data$country.models)[i],
                               "p" = data$country.models[[i]]$specs$lags$lags$domestic,
                               "p.star" = data$country.models[[i]]$specs$lags$lags$foreign,
                               "s.global" = data$country.models[[i]]$specs$lags$lags$global,
                               "rank" = data$country.models[[i]]$specs$rank$rank)
                    )
  }
  return(result)
}