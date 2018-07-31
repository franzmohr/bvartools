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
  if (any(class(data) == "list")) {
    if(sum(unlist(lapply(data, class)) == "ts") != length(data)) {stop("Data must be of class 'ts'.")} 
    
    f <- function(x, variables){
      if (is.null(variables)){
        for (i in dimnames(x)[[2]]){
          x[, i] <- c(NA, diff(x[, i]))
        }
      } else {
        for (i in variables){
          if (is.element(i, dimnames(x)[[2]])){
            x[, i] <- c(NA, diff(x[, i]))
          }
        } 
      }
      x <- stats::na.omit(x)
      attributes(x)$na.action <- NULL
      return(x)
    }
    
    data <- lapply(data, f, variables)
  } else {
    if(!"ts" %in% class(data)) {stop("Data must be of class 'ts'.")}
    
    if (is.null(variables)){
      for (i in dimnames(data)[[2]]){
        data[, i] <- c(NA, diff(data[, i]))
      }
    } else {
      for (i in variables){
        if (is.element(i, dimnames(data)[[2]])){
          data[, i] <- c(NA, diff(data[, i]))
        }
      } 
    }
    data <- stats::na.omit(data)
    attributes(data)$na.action <- NULL
  }
  return(data)
}