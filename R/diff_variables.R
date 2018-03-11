#' Differences of Variables
#' 
#' Produces the first difference of selected time series in a list of data.
#' 
#' @param data a "zoo" object or a list of named "zoo" objects.
#' @param variables a character vector of variables that should be differenced. If \code{NULL} (default), all variables
#' will be differenced.
#' 
#' @return A differenced "zoo" object or a named list of differenced "zoo" objects.
#' 
#' @examples
#' \dontrun{
#' data(gvar2013)
#' 
#' country.data <- gvar2013$country.data
#' 
#' # Only variable y and q are differenced
#' country.data.yq <- diff_variables(data = country.data, variables = c("y", "q"))
#' 
#' # All variables are differenced
#' country.data.all <- diff_variables(data = country.data)
#' }
#' 
#' @export
diff_variables <- function(data, variables = NULL){
  if(!requireNamespace("zoo")) {stop("Function requires the package 'zoo'.")}
  
  if (class(data) == "list") {
    if(any(unlist(lapply(data, class)) != "zoo")) {stop("Data must be of class 'zoo'.")} 
    
    f <- function(x, variables){
      if (is.null(variables)){
        for (i in names(x)){
          x[, i] <- diff(x[, i], na.pad = TRUE)
        }
      } else {
        for (i in variables){
          if (is.element(i, names(x))){
            x[, i] <- diff(x[, i], na.pad = TRUE)
          }
        } 
      }
      x <- stats::na.omit(x)
      attributes(x)$na.action <- NULL
      return(x)
    }
    
    data <- lapply(data, f, variables)  
  } else {
    if(class(data) != "zoo") {stop("Data must be of class 'zoo'.")}
    
    if (is.null(variables)){
      for (i in names(data)){
        data[, i] <- diff(data[, i], na.pad = TRUE)
      }
    } else {
      for (i in variables){
        if (is.element(i, names(data))){
          data[, i] <- diff(data[, i], na.pad = TRUE)
        }
      } 
    }
    data <- stats::na.omit(data)
    attributes(data)$na.action <- NULL
  }
  
  return(data)
}