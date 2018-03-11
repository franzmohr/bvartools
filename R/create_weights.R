#' Create a Weight Matrix
#' 
#' Produces a matrix of fixed weights or an array for time varying weights for a GVAR model.
#' 
#' @param weight.data a matrix or an array containing weight data used to generate country-specific weight matrices.
#' @param period a character vector or a numeric specifying the periods used to calculate the final weights, see details below.
#' @param country.data a named list of "zoo" objects containing country data. Only used if \code{period} is a numeric.
#' 
#' @details
#' If a character vector is provided for \code{period}, the function will calculate country weights based on the sums over the
#' specified periods. If a numeric is proved, the weights will be construced from rolling sums over the last \code{period} periods.
#' If the country data starts earlier than the weight series, the sums over the first \code{period} observations
#' of the weight data will be used until the periods match.
#' 
#' @return A matrix or an array containing weight data used to generate country-specific weight matrices. Each row sums up to 1.
#' 
#' @examples
#' \dontrun{
#' # Constant weight matrix 
#' data(gvar2016)
#' 
#' weight.data <- gvar2016$weight.data
#' period <- c("2008", "2009", "2010", "2011", "2012")
#' 
#' weight.data <- create_weights(weight.data, period)
#' 
#' 
#' # Time varying weight matrix
#' data(gvar2016)
#' 
#' weight.data <- gvar2016$weight.data
#' country.data <- gvar2016$country.data
#' period <- 3
#' 
#' weight.data <- create_weights(weight.data, period, country.data)
#' }
#' @export
create_weights <- function(weight.data, period, country.data = NULL){
  if (is.null(dimnames(weight.data)[[1]]) | is.null(dimnames(weight.data)[[2]])) {
    stop("Weight matrix rows and columns must both be named.")
  }
  if (any(apply(weight.data, 3, diag) != 0)) {stop("The diagonal elements of the weight matrix must be zero.")}
  if (class(period) == "character") {
    w <- weight.data[,,period]
    w <- apply(w, c(1, 2), sum)
    w <- t(apply(w, 1, function(x){return(x/sum(x))}))
  }
  if (class(period) == "numeric") {
    if (is.null(country.data)){stop("The construction of time varying weights requries to specify the argument country.data.")}
    tt <- unique(unlist(lapply(country.data, NROW)))
    if (length(tt) > 1) {stop("Objects in country.data do not have the same number of observations.")}
    t.ind <- zoo::index(country.data[[1]])
    t.avail <- as.numeric(dimnames(weight.data)[[3]])
    if (class(t.ind) == "yearqtr"){
      t.temp <- as.numeric(substring(as.character(t.ind), 1, 4))
    } else {
      stop("Sorry, the 'bgvars' package currently only supports the 'yearqtr' time series format.")
    }
    t <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (t.temp[i] <= t.avail[period]) {
        t[i,] <- t.avail[1:period]
      }
      if (t.temp[i]%in%t.avail[(period + 1):length(t.avail)]) {
        t[i,] <- (t.temp[i] - period + 1):t.temp[i]
      }
    }
    w <- array(NA, dim = c(dim(weight.data)[1], dim(weight.data)[2], tt))
    dimnames(w)[[1]] <- dimnames(weight.data)[[1]]
    dimnames(w)[[2]] <- dimnames(weight.data)[[2]]
    for (i in 1:tt){
      w.temp <- apply(weight.data[,,as.character(t[i,])], c(1, 2), sum)
      w[,,i] <- t(apply(w.temp, 1, function(x) {return(x/sum(x))}))
    }
    dimnames(w)[[3]] <- t.ind
  }
  return(w)
}
