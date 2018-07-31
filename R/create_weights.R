#' Create a Weight Matrix
#' 
#' Produces a matrix of fixed weights or an array of time varying weights from raw data for a GVAR model.
#' 
#' @param weight.data a list of named time series or a named matrix containing data for the construction of weights.
#' @param period an integer for time varying weights or a numeric vector specifying the periods used to calculate constant weights (see 'Details').
#' @param country.data a named list of time series containing country data. Only requried if \code{period} is an integer.
#' 
#' @details
#' The function assists in the creation of country-specific weight matrices from raw data. If a numeric vector
#' is provided for \code{period}, the function will calculate country weights based on the sums over the
#' specified periods. If an integer is proved, the weights will be construced from rolling sums over the last \code{period} periods.
#' If the country data starts earlier than the weight series, the sums over the first \code{period} observations
#' of \code{weight.data} will be used until the periods match.
#' 
#' @return A named matrix or array containing country-specific weight matrices.
#' 
#' @examples
#' # Constant weights 
#' data(gvar2016)
#' 
#' weight.data <- gvar2016$weight.data
#' period <- 2008:2012
#' 
#' weight.data <- create_weights(weight.data, period)
#' 
#' 
#' # Time varying weights
#' data(gvar2016)
#' 
#' weight.data <- gvar2016$weight.data
#' country.data <- gvar2016$country.data
#' period <- 3
#' 
#' weight.data <- create_weights(weight.data, period, country.data)
#' 
#' @export
create_weights <- function(weight.data, period, country.data = NULL){
  # Check if weight.data is a list or matrix.
  if (class(weight.data) %in% c("matrix", "list")) {
    if (class(weight.data) == "list") {
      l <- TRUE
    } else {
      l <- FALSE
    } 
  } else {
    stop("weight.data must be of class list or matrix.")
  }
  
  if (!l) {
    if (is.null(dimnames(weight.data)[[1]]) | is.null(dimnames(weight.data)[[2]])) {
      stop("Both the rows and columns of the weight.matrix must be named.")
    }
    if (any(apply(weight.data, 3, diag) != 0)) {stop("The diagonal elements of the weight.matrix must be zero.")}
  }
  
  if (length(period) > 1) {
    period <- as.character(period)
    if (l){
      w <- matrix(NA, length(weight.data), length(weight.data))
      dimnames(w) <- list(names(weight.data), names(weight.data))
      for (i in names(weight.data)) {
        temp <- colSums(weight.data[[i]][period, ])
        w[i, names(temp)] <- temp / sum(temp)
      }
    }
  } else {
    if (is.null(country.data)){stop("The construction of time varying weights requries to specify the argument country.data.")}
    tt <- unique(unlist(lapply(country.data, NROW)))
    if (length(tt) > 1) {stop("Objects in country.data do not have the same number of observations.")}
    t.temp <- as.numeric(time(country.data[[1]]))
    
    t.ind <- tsp(country.data[[1]])
    t.avail <- unique(unlist(lapply(weight.data, NROW)))
    if (length(t.avail) > 1) {stop("Objects in weight.data do not have the same number of observations.")}
    t.avail <- as.numeric(time(weight.data[[1]]))

    t <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (t.temp[i] <= t.avail[period]) {
        t[i,] <- t.avail[1:period]
      }
      if (t.temp[i] >= t.avail[period]) {
        pos_t <- which(floor(t.temp[i]) == t.avail)
        pos_t <- (pos_t - period + 1):pos_t
        t[i,] <- t.avail[pos_t]
      }
    }
    
    w <- array(NA, dim = c(length(weight.data), length(weight.data), tt))
    dimnames(w) <- list(names(weight.data), names(weight.data), as.character(t.temp))
    for (i in 1:tt){
      for (j in names(country.data)) {
        temp <- colSums(weight.data[[j]][rownames(weight.data[[j]]) %in% as.character(t[i,]), ])
        w[j, names(temp), i] <- temp / sum(temp)
      }
    }
  }
  return(w)
}
