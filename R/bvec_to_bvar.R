#' Transform a VECM to VAR in levels
#' 
#' An object of class `bvec`` is transformed to a VAR in level presentation.
#' 
#' @param object an object of class `bvec`.
#' 
#' @return An object of class `bvar`.
#' 
#' @export
bvec_to_bvar <- function(object) {
  draws <- nrow(object$Pi)
  k <- NROW(object$y)
  
  if (!is.null(object$Gamma)) {
    p <- NCOL(object$Gamma) / k^2
    p <- p + 1
    W <- diag(-1, k * p)
    W[1:k, 1:k] <- diag(1, k)
    W[-(1:k), -(k * (p - 1) + 1:k)] <- W[-(1:k),-(k * (p - 1) + 1:k)] + diag(k * (p - 1))
    J <- matrix(0, k, k * p)
    J[1:k, 1:k] <- diag(1, k)
    
    A <- matrix(NA, k^2 * p, draws)
    for (draw in 1:draws) {
      A[, draw] <- cbind(matrix(object$Pi[draw, ], k), matrix(object$Gamma[draw, ], k)) %*% W + J
    }
  } else {
    A <- matrix(NA, k^2, draws)
    for (draw in 1:draws) {
      A[, draw] <- matrix(object$Pi[draw, ], k) + matrix(diag(1, k), k)
    }
  }

  if (!is.null(object$Ypsilon)) {
    k_x <- nrow(object$ect) - k
    s <- NCOL(object$Ypsilon) / (k * k_x)
    s <- s + 1
    W <- diag(-1, k_x * (s + 1))
    W[1:k_x, 1:k_x] <- 0
    W[1:k_x, k_x + 1:k_x] <- diag(1, k_x)
    W[-(1:k_x), 1:(k_x * s)] <- W[-(1:k_x), 1:(k_x * s)] + diag(1, k_x * s)
    
    B <- matrix(NA, k * k_x * (s + 1), draws)
    for (draw in 1:draws){
        B[, draw] <- cbind(matrix(object$Pi_x[draw, ], k), matrix(object$Ypsilon[draw, ], k)) %*% W
    }
  } else {
    B <- NULL
  }
  
  if (!is.null(object$y)) {
    y <- object$y
  } else {
    y <- NULL
  }
  if (!is.null(object$x)) {
    x <- object$x
  } else {
    x <- NULL
  }
  if (!is.null(object$D)) {
    D <- t(object$D)
  } else {
    D <- NULL
  }
  if (!is.null(object$Sigma)) {
    Sigma <- t(object$Sigma)
  } else {
    Sigma <- NULL
  }
  
  dimnames(y)[[1]] <- gsub("d.", "", dimnames(y)[[1]])
  
  result <- bvar(y = y, x = x, A = A, B = B, D = D, Sigma = Sigma)
  return(result)
}