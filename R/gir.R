#' Generalised Impulse Response
#' 
#' Produces the generalised impulse responses of a moving average representation of a GVAR model.
#' 
#' @param G matrix or array containing the coefficients of a VAR model.
#' @param G0.i Inverse of G0 from the GVAR model.
#' @param Sigma Variance-covariance matrix of the GVAR model.
#' @param n.ahead an integer specifying the steps.
#' @param impulse Character vector with two elements specifying the impulse country and variable.
#' @param response Character vector with two elements specifying the response country and variable.
#' @param index Data frame containing indices of countries and variables in the model.
#' @param shock Size of shock. Either "sd" (default), "nsd" for negative standard deviation of the impulse variable or a numeric value.
#' 
#' @return Vector of generalised impulse response.
#' 
#' @export
gir <- function(G, G0.i, Sigma, n.ahead, impulse, response, index, shock = "sd", t = NULL){
  k <- dim(G)[1]
  p <- dim(G)[2] / k
  if (n.ahead < p) {p <- n.ahead}
  
  par <- c()
  for (i in 1:n.ahead){
    par <- c(par, list(matrix(0, k, k)))
  }
  for (i in 1:p){
    par[[i]] <- G[, (i - 1) * k + 1:k]
  }
  
  # Generate MA representation
  Phi <- list(diag(1, k))
  for (i in 1:n.ahead){
    Phi <- c(Phi, list(matrix(0, k, k)))
  }
  for (i in 1:n.ahead){
    Phi.i <- matrix(0, k, k)
    for (j in 1:i){
      Phi.i <- Phi.i + Phi[[i - j + 1]]%*%par[[j]]
    }
    Phi[[i + 1]] <- Phi.i
  }
  
  impulse <- which(index[, 1] == impulse[1] & index[, 2] == impulse[2])
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(index[, 1] == response[1] & index[, 2] == response[2])
  if (length(response) == 0){stop("Response variable not available.")}
  
  sjj <- Sigma[impulse, impulse]
  sjj.sqrt <- sqrt(sjj)
  if (shock == "sd") {
    s <- 1
  } else {
    if (shock == "nsd"){
      s <- -1
    } else {
      s <- shock / sjj.sqrt
    }
  }
  ej <- matrix(0, k, 1)
  ej[impulse, 1] <- 1
  GI <- unlist(lapply(Phi,
                      function(x, G0.i, Sigma, ej, sjj.sqrt, s, response) {
                        return(((x%*%G0.i%*%Sigma%*%ej)/sjj.sqrt*s)[response, 1])
                      }, G0.i, Sigma, ej, sjj.sqrt, s, response)
  )
  return(GI)
}