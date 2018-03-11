#' Obtain Country Model Data
#' 
#' Multiplies the country-specific weight matrix with the data of the global model to obtain domestic
#' and foreign variables of a country model.
#' 
#' @param W a country model's weight matrix as produced by \code{\link{weight_matrices}}.
#' @param X a "zoo" object containing the data of the global model.
#' 
#' @return A matrix.
#' 
#' @export
get_z <- function(W, X){
  X <- t(X)
  if (class(W) == "matrix"){
    x <- W%*%X
  } else {
    x <- matrix(NA, dim(W)[1], dim(W)[3])
    for (i in 1:dim(W)[3]){
      x[, i] <- W[,, i]%*%X[, i]
    } 
  }
  dimnames(x) <- dimnames(W)[1]
  return(x)
}