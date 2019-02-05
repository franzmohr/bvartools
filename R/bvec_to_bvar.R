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
    m <- NCOL(object$Pi_x)
    s <- NCOL(object$Ypsilon) / (k * m)
    W <- diag(-1, m * (s + 1))
    W[1:m, 1:m] <- 0
    W[1:m, m + 1:m] <- diag(1, m)
    W[-(1:m), 1:(m * s)] <- W[-(1:m), 1:(m * s)] + diag(1, m * s)
    
    B <- matrix(NA, k * m * (s + 1), draws)
    for (draw in 1:draws){
        B[, draw] <- cbind(matrix(object$Pi_x[draw, ], k), matrix(object$Ypsilon[draw, ], k)) %*% W
    }
  } else {
    B <- NULL
    m <- 0
    s <- 0
  }
  
  y <- object$y
  y_names <- gsub("d.", "", dimnames(y)[[1]])
  dimnames(y) <- list(y_names, NULL)
  if (!is.null(object$w)) {
    y <- y + matrix(object$w[1:k, ], k)
    dimnames(y) <- list(y_names, NULL)
  }
  
  if (!is.null(object$x)) {
    if (p > 1) {
      y_diff <- matrix(object$x[1:(k * (p - 1)), 1], k)
      y_diff <- matrix(y_diff[,(p - 1):1], k)
      y_diff <- cbind(y_diff, object$y[, 1])
      y_level <- y_diff * NA
      y_level[, p] <- y[, 1] - y_diff[, p]
      for (i in (p - 1):1){
        y_level[, i] <- y_level[, i + 1] - y_diff[, i]
      }
      y <- stats::ts(t(cbind(y_level, y)))
    } else {
      if (!is.null(object$w)) {
        y <- stats::ts(t(cbind(matrix(object$w[1:k, 1], k), y)))
      }
    }
    
    if (m > 0) {
      x_names <- dimnames(object$x)[[1]][k * (p - 1) + 1:m]
      x_names <- gsub("d.", "", x_names)
      x_names <- gsub(".0", "", x_names)
      if (s > 0) {
        x0 <- matrix(object$x[(k * (p - 1)) + 1:m,], m)
        x0 <- x0 + matrix(object$w[(k * (p - 1)) + 1:m, ], m)
        x_diff <- matrix(object$x[(k * (p - 1)) + 1:(m * s), 1], m)
        x_diff <- matrix(x_diff[, s:1], m)
        x_level <- x_diff * NA
        x_level[, s] <- x0[, 1] - x_diff[, s]
        for (i in (s - 1):1){
          x_level[, i] <- x_level[, i + 1] - x_diff[, i]
        }
        x <- stats::ts(t(cbind(x_level, x0)))
      }
    }
    
    temp <- y
    temp_names <- y_names
    for (i in 1:p) {
      temp <- cbind(temp, stats::lag(y, -i))
      temp_names <- c(temp_names, paste(y_names, i, sep = "."))
    }
    if (m > 0) {
      temp <- cbind(temp, x)
      temp_names <- c(temp_names, paste(x_names, 0, sep = "."))
      for (i in 1:s) {
        temp <- cbind(temp, stats::lag(x, -i))
        temp_names <- c(temp_names, paste(x_names, i, sep = "."))
      } 
    }
    
    temp <- stats::na.omit(temp)
    dimnames(temp)[[2]] <- temp_names
    y <- t(temp[, 1:k])
    dimnames(y) <- list(y_names, NULL)
    x <- t(temp[, -(1:k)])
    dimnames(x) <- list(temp_names[-(1:k)], NULL)
  } else {
    x <- NULL
  }
  
  if (!is.null(object$A0)) {
    A0 <- t(object$A0)
  } else {
    A0 <- NULL
  }
  
  n <- 0
  n_Pi_d <- 0
  C <- NULL
  x_det <- NULL
  x_det_names <- NULL
  if (!is.null(object$Pi_d)) {
    n <- NCOL(object$Pi_d) / k
    n_Pi_d <- n
    C <- t(object$Pi_d)
    if (!is.null(object$w)) {
      x_det_names <- dimnames(object$w)[[1]][-(1:(k + m))]
      x_det <- matrix(object$w[-(1:(k + m)),], n)
    }
  }
  
  n_c <- 0
  if (!is.null(object$C)) {
    n_c <- NCOL(object$C) / k
    n <- n + n_c
    if (!is.null(C)) {
      C <- rbind(C, t(object$C))
    } else {
      C <- t(object$C)
    }
  }
  
  if (!is.null(object$x)) {
    x_det_names <- c(x_det_names, dimnames(object$x)[[1]][-(1:(k * (p - 1) + m * s))])
    x_temp <- matrix(object$x[-(1:(k * (p - 1) + m * s)),], n_c)
    x_det <- rbind(x_det, x_temp)
    dimnames(x_det) <- list(x_det_names, NULL)
    x <- rbind(x, x_det)
  }
  
  if (!is.null(object$Sigma)) {
    Sigma <- t(object$Sigma)
  } else {
    Sigma <- NULL
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
  
  result <- bvar(data = data, exogen = exogen, y = y, x = x, A0 = A0, A = A, B = B, C = C, Sigma = Sigma)
  return(result)
}