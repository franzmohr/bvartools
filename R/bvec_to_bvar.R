#' Transform a VECM to VAR in levels
#' 
#' An object of class "bvec" is transformed to a VAR in level representation.
#' 
#' @param object an object of class "bvec".
#' 
#' @return An object of class "bvar".
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
    k_x <- NCOL(object$Pi_x)
    s <- NCOL(object$Ypsilon) / (k * k_x)
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
  
  if (!is.null(object$data)) {
    data <- object$data
  } else {
    data <- NULL
  }
  
  if (!is.null(object$exogen)) {
    exogen <- object$exogen
  } else {
    exogen <- NULL
  }
  
  if (!is.null(object$y)) {
    y <- object$y
    dimnames(y)[[1]] <- gsub("d.", "", dimnames(y)[[1]])
  } else {
    y <- NULL
  }
  if (!is.null(object$x)) {
    x <- object$x
    dimnames(x)[[1]] <- gsub("d.", "", dimnames(x)[[1]])
  } else {
    x <- NULL
  }
  
  if (!is.null(object$A0)) {
    A0 <- t(object$A0)
  } else {
    A0 <- NULL
  }
  
  if (!is.null(object$Pi_d)) {
    Pi_d <- t(object$Pi_d)
  } else {
    Pi_d <- NULL
  }
  if (!is.null(object$C)) {
    if (!is.null(Pi_d)) {
      C <- rbind(Pi_d, t(object$C))
    } else {
      C <- t(object$C)
    }
  } else {
    C <- NULL
    if (!is.null(Pi_d)) {
      C <- Pi_d
    }
  }
  if (!is.null(object$Sigma)) {
    Sigma <- t(object$Sigma)
  } else {
    Sigma <- NULL
  }
  
  result <- bvar(data = data, exogen = exogen, y = y, x = x, A0 = A0, A = A, B = B, C = C, Sigma = Sigma)
  return(result)
}