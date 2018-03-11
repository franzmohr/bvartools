#' Obtain Global Model Data
#' 
#' Produces data for the global VAR model and updates country model specifications.
#' 
#' @param country.data a named list of "zoo" objects containing the country data.
#' @param country.specs a list of model specifications as produced by \code{\link{country_specifications}}.
#' @param global.data a "zoo" object with global data.
#'
#' @details The function extracts the specified variables of each country model from the original data source
#' and adds them to the time series of the endogenous variables of global model. If a country model 
#' specification includes a common global variable as an endogenous domestic variable, it will be added to the endogenous global
#' sample and dropped from the time series of common global variables.
#'
#' @return A list containing the final global time series of endogenous variables, exogenous global variables,
#' an index of variables and final country specifications:
#' \item{endogenous.variables}{a "zoo" object containing the endogenous variables of the global model.}
#' \item{global.data}{a "zoo" object containing common global variables.}
#' \item{index}{a data frame with country and variable names.}
#' \item{country.specs}{a list of updated country model specifications.}
#' 
#' @export
global_series <- function(country.data, country.specs, global.data = NULL){
  n.c <- length(country.specs) # Number of countries
  c.names <- names(country.specs) # Country names
  # Get list of all variables
  avail.dom.vars <- unique(unlist(lapply(country.data, names))) # Domestic variables
  if (!is.null(global.data)){
    avail.global.vars <- names(global.data) # Global variables
  }
  
  # Generate final variable index
  ind <- data.frame() # Index
  X <- zoo::zoo() # Final data set of domestic variables (incl. global endogenous)
  time.class <- unique(unlist(lapply(country.data, function(x) {return(class(zoo::index(x)))})))
  if (length(time.class) > 1) {stop("Country data must have the same class.")}
  class(zoo::index(X)) <- time.class
  X.names <- c()
  for (i in c.names){
    temp <- zoo::zoo()
    class(zoo::index(temp)) <- time.class
    temp.names <- c()
    for (j in country.specs[[i]]$domestic.variables){
      if (is.element(j, names(country.data[[i]]))){
        temp <- cbind(temp, country.data[[i]][, j])
      } else {
        temp <- cbind(temp, global.data[, j])
      }
      temp.names <- append(temp.names, j)
    }
    X <- cbind(X, temp)
    n.temp <- length(temp.names)
    X.names <- c(X.names, temp.names)
    ind <- rbind(ind, data.frame("Country" = rep(i, n.temp),
                                 "Variable" = temp.names))
  }
  X <- stats::na.omit(X)
  names(X) <- X.names
  ind[, 1] <- as.character(ind[, 1])
  ind[, 2] <- as.character(ind[, 2])
  
  # Update country specifications to correct for endogenous global variables
  variables <- unique(ind[, 2]) # All endogenous variables
  endo.g.var <- NULL
  for (i in c.names){
    # Check if there is an endogenous global variable in a country specification
    if (any(is.element(country.specs[[i]]$global.variables, variables))){
      star <- country.specs[[i]]$star.variables
      global <- country.specs[[i]]$global.variables
      pos.endo.g <- which(is.element(global, variables))
      endo.global <- global[pos.endo.g]
      endo.g.var <- c(endo.g.var, endo.global)
      star <- c(star, endo.global)
      global <- global[-pos.endo.g]
      if (length(global) == 0){
        global <- NA
        country.specs[[i]]$lags$lags$global <- NA
      }
      country.specs[[i]]$star.variables <- star
      country.specs[[i]]$global.variables <- global
    }
  }
  
  # Update global data
  endo.g.var <- unique(endo.g.var)
  if (length(endo.g.var) > 0){
    if (length(endo.g.var) == length(names(global.data))){
      global.data <- NULL
    } else {
      global.data <- global.data[, -which(is.element(endo.g.var, names(global.data)))]
    }
  }
  
  result <- list(endogenous.variables = X, global.data = global.data,
                 index = ind, country.specs = country.specs) 
  return(result)
}