#' Orthogonalised and Structural Impulse Responses
#' 
#' The function produces a draw of an orthogonalised or structural impulse response.
#' 
#' @param x a list.
#' @param impulse integer specifying the impulse variable.
#' @param response integer specifying the response variable.
#' @param type character specifying the type of shock.
#' 
#' @return A numeric vector.
#' 
#' @export
sir <- function(x, impulse, response, type) {
  Phi <- x$Phi
  if (type == "oir") {
    P <- t(chol(x$Sigma))
  }
  if (type == "sir") {
    P <- x$A0
  }
  
  n <- NCOL(Phi)
  h <- nrow(Phi) / n
  
  result <- rep(NA, h)
  for (i in 1:h) {
    result[i] <- (Phi[(i - 1) * n + 1:n,] %*% P)[response, impulse]
  }
  return(result)
}