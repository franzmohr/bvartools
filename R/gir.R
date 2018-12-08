#' Generalised Impulse Responses
#' 
#' The function produces a draw of an orthogonalised or structural impulse response.
#' 
#' @param x a list.
#' @param impulse integer specifying the impulse variable.
#' @param response integer specifying the response variable.
#' @param shock character or numeric specifying the type of shock.
#' 
#' @return A numeric vector.
#' 
#' @references
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters}, 58, 17-29.
#' 
#' @export
gir <- function(x, impulse, response, shock) {
  Phi <- x$Phi
  Sigma <- x$Sigma
  sjj_sqrt <- sqrt(Sigma[impulse, impulse])
  
  n <- NROW(Sigma)
  h <- nrow(Phi) / n
  
  if (shock == "sd") {
    s <- 1
  } else {
    if (shock == "nsd"){
      s <- -1
    } else {
      s <- shock / sjj_sqrt
    }
  }
  
  ej <- matrix(0, n, 1)
  ej[impulse, 1] <- 1
  
  result <- rep(NA, h)
  for (i in 1:h) {
    result[i] <- (Phi[(i - 1) * n + 1:n,] %*% Sigma %*% ej / sjj_sqrt * s)[response, 1]
  }
  return(result)
}