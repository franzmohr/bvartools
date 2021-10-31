#' Transform a VEC Model to a VAR in Levels
#' 
#' An object of class \code{"bvec"} is transformed to a VAR in level representation.
#' 
#' @param object an object of class \code{"bvec"}.
#' 
#' @return An object of class \code{"bvar"}.
#' 
#' @examples
#' 
#' # Load data
#' data("e6")
#' # Generate model
#' data <- gen_vec(e6, p = 4, r = 1, const = "unrestricted", season = "unrestricted")
#' # Obtain data matrices
#' y <- t(data$data$Y)
#' w <- t(data$data$W)
#' x <- t(data$data$X)
#' 
#' # Reset random number generator for reproducibility
#' set.seed(1234567)
#' 
#' iterations <- 400 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations should be much higher, e.g. 30000.
#' 
#' burnin <- 100 # Number of burn-in draws
#' draws <- iterations + burnin
#' 
#' r <- 1 # Set rank
#' 
#' tt <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' k_w <- nrow(w) # Number of regressors in error correction term
#' k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
#' 
#' k_alpha <- k * r # Number of elements in alpha
#' k_beta <- k_w * r # Number of elements in beta
#' k_gamma <- k * k_x
#' 
#' # Set uninformative priors
#' a_mu_prior <- matrix(0, k_x * k) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, k_x * k) # Inverse of the prior covariance matrix
#' 
#' v_i <- 0
#' p_tau_i <- diag(1, k_w)
#' 
#' u_sigma_df_prior <- r # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
#' u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
#' 
#' # Initial values
#' beta <- matrix(c(1, -4), k_w, r)
#' u_sigma_i <- diag(1 / .0001, k)
#' g_i <- u_sigma_i
#' 
#' # Data containers
#' draws_alpha <- matrix(NA, k_alpha, iterations)
#' draws_beta <- matrix(NA, k_beta, iterations)
#' draws_pi <- matrix(NA, k * k_w, iterations)
#' draws_gamma <- matrix(NA, k_gamma, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#' 
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'   # Draw conditional mean parameters
#'   temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = u_sigma_i,
#'                          v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
#'                          gamma_mu_prior = a_mu_prior,
#'                          gamma_v_i_prior = a_v_i_prior)
#'   alpha <- temp$alpha
#'   beta <- temp$beta
#'   Pi <- temp$Pi
#'   gamma <- temp$Gamma
#'   
#'   # Draw variance-covariance matrix
#'   u <- y - Pi %*% w - matrix(gamma, k) %*% x
#'   u_sigma_scale_post <- solve(tcrossprod(u) +
#'      v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
#'   u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'   u_sigma <- solve(u_sigma_i)
#'   
#'   # Update g_i
#'   g_i <- u_sigma_i
#'   
#'   # Store draws
#'   if (draw > burnin) {
#'     draws_alpha[, draw - burnin] <- alpha
#'     draws_beta[, draw - burnin] <- beta
#'     draws_pi[, draw - burnin] <- Pi
#'     draws_gamma[, draw - burnin] <- gamma
#'     draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#' 
#' # Number of non-deterministic coefficients
#' k_nondet <- (k_x - 4) * k
#' 
#' # Generate bvec object
#' bvec_est <- bvec(y = data$data$Y, w = data$data$W,
#'                  x = data$data$X[, 1:6],
#'                  x_d = data$data$X[, 7:10],
#'                  Pi = draws_pi,
#'                  Gamma = draws_gamma[1:k_nondet,],
#'                  C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
#'                  Sigma = draws_sigma)
#' 
#' # Thin posterior draws
#' bvec_est <- thin_posterior(bvec_est, thin = 5)
#' 
#' # Transfrom VEC output to VAR output
#' bvar_form <- bvec_to_bvar(bvec_est)
#' 
#' 
#' @export
bvec_to_bvar <- function(object) {
  
  if (!any(class(object) %in% "bvec")) {
    stop("Argument 'object' must be of class 'bvec'.")
  }
  
  draws <- NULL
  specs <- NULL
  vars <- c("Pi", "Pi_x", "Pi_d", "Gamma", "Upsilon", "C", "A0")
  for (i in vars) {
    if (is.null(draws)) {
      if (!is.null(object[[i]])) {
        if (is.list(object[[i]])) {
          draws <- nrow(object[[i]][[1]])
        } else {
          draws <- nrow(object[[i]]) 
        }
      }
    }
    if (is.null(specs)) {
      if (is.list(object[[i]])) {
        specs <- attr(object[[i]][[1]], "mcpar")
      } else {
        specs <- attr(object[[i]], "mcpar")
      }
    }
  }
  
  k <- NCOL(object[["y"]])
  tt <- NROW(object[["y"]])
  tvp <- object[["specifications"]][["tvp"]]
  p <- object[["specifications"]][["lags"]][["p"]]
  r <- object[["specifications"]][["rank"]]
  
  A <- NULL
  if (!is.null(object[["Gamma"]])) {
    
    W <- diag(-1, k * p)
    W[1:k, 1:k] <- diag(1, k)
    W[-(1:k), -(k * (p - 1) + 1:k)] <- W[-(1:k),-(k * (p - 1) + 1:k)] + diag(k * (p - 1))
    J <- matrix(0, k, k * p)
    J[1:k, 1:k] <- diag(1, k)
    
    n_gamma <- k * k * p
    
    if (!tvp[["Pi"]] & !tvp[["Gamma"]]) {
      
      A <- matrix(NA, n_gamma, draws)
      for (draw in 1:draws) {
        if (is.null(object[["Pi"]])) {
          A[, draw] <- cbind(matrix(0, k, k), matrix(object[["Gamma"]][draw, ], k)) %*% W + J
        } else {
          A[, draw] <- cbind(matrix(object[["Pi"]][draw, ], k), matrix(object[["Gamma"]][draw, ], k)) %*% W + J 
        }
      } 
      
    } else {
      
      A <- matrix(NA, n_gamma * tt, draws)
      pi_temp <- matrix(0, k, k)
      for (draw in 1:draws) {
        for (i in 1:tt) {
          if (!is.null(object[["Pi"]])) {
            if (tvp[["Pi"]]) {
              pi_temp <- matrix(object[["Pi"]][[i]][draw, ], k)
            } else {
              pi_temp <- matrix(object[["Pi"]][draw, ], k)
            }
          }
          if (tvp[["Gamma"]]) {
            gamma_temp <- matrix(object[["Gamma"]][[i]][draw, ], k)
          } else {
            gamma_temp <- matrix(object[["Gamma"]][draw, ], k)
          }
          A[(i - 1) * n_gamma + 1:n_gamma, draw] <- cbind(pi_temp, gamma_temp) %*% W + J 
        }
      }
      
    }
    
  } else {
    if (!is.null(object[["Pi"]])) {
      n_a <- k * k
      if (tvp[["Pi"]]) {
        A <- matrix(NA, n_a * tt, draws)
        for (draw in 1:draws) {
          for (i in 1:tt) {
            A[(i - 1) * n_a + 1:n_a, draw] <- matrix(object[["Pi"]][[i]][draw, ], k) + matrix(diag(1, k), k)
          }
        }
      } else {
        A <- matrix(NA, n_a, draws)
        for (draw in 1:draws) {
          A[, draw] <- matrix(object[["Pi"]][draw, ], k) + matrix(diag(1, k), k)
        } 
      }
    } else {
      A <- matrix(0, k * k, draws)
    }
  }
  
  B <- NULL
  m <- 0
  s <- 0
  if (!is.null(object[["Upsilon"]])) {
    
    m <- object[["specifications"]][["dims"]][["M"]]
    s <- object[["specifications"]][["lags"]][["s"]]
    
    W <- diag(-1, m * (s + 1))
    W[1:m, 1:m] <- 0
    W[1:m, m + 1:m] <- diag(1, m)
    W[-(1:m), 1:(m * s)] <- W[-(1:m), 1:(m * s)] + diag(1, m * s)
    
    n_b <- k * m * (s + 1)
    
    if (!tvp[["Pi_x"]] & !tvp[["Upsilon"]]) {
      
      B <- matrix(NA, n_b, draws)
      for (draw in 1:draws){
        if (!is.null(object[["Pi_x"]])) {
          B[, draw] <- cbind(matrix(object[["Pi_x"]][draw, ], k), matrix(object[["Upsilon"]][draw, ], k)) %*% W 
        } else {
          B[, draw] <- cbind(matrix(0, k, m), matrix(object[["Upsilon"]][draw, ], k)) %*% W 
        }
      } 
      
    } else {
      
      B <- matrix(NA, n_b * tt, draws)
      pix_temp <- matrix(0, k, m)
      for (draw in 1:draws){
        for (i in 1:tt) {
          if (!is.null(object[["Pi_x"]])) {
            if (tvp[["Pi_x"]]) {
              pix_temp <- matrix(object[["Pi_x"]][[i]][draw, ], k)
            } else {
              pix_temp <- matrix(object[["Pi_x"]][draw, ], k)
            } 
          }
          if (tvp[["Upsilon"]]) {
            ups_temp <- matrix(object[["Upsilon"]][[i]][draw, ], k)
          } else {
            ups_temp <- matrix(object[["Upsilon"]][draw, ], k)
          }
          B[(i - 1) * n_b + 1:n_b, draw] <- cbind(pix_temp, ups_temp) %*% W 
        }
      }
    }
  }
  
  y <- object[["y"]]
  y_names <- gsub("d.", "", dimnames(y)[[2]])
  dimnames(y) <- list(NULL, y_names)
  # Reconstruct y if w is available
  if (!is.null(object[["w"]])) {
    y <- y + as.matrix(object[["w"]][, 1:k])
    dimnames(y) <- list(NULL, y_names)
  }
  tsp_temp <- stats::tsp(y)
  
  if (p > 1) {
    y_diff <- matrix(object[["x"]][1,], k) # Get first observations
    y_diff <- matrix(y_diff[,(p - 1):1], k) # Reverse order
    y_diff <- cbind(y_diff, matrix(object[["y"]][1, ], k))
    y_level <- y_diff * NA
    y_level[, p] <- y[1, ] - y_diff[, p]
    for (i in (p - 1):1){
      y_level[, i] <- y_level[, i + 1] - y_diff[, i]
    }
    y <- rbind(t(y_level), y)
  } else {
    if (!is.null(object[["w"]])) {
      y <- rbind(object[["w"]][1, ], y)
    }
  }
  
  if (m > 0) {
    x_names <- dimnames(object[["x_x"]])[[2]][1:m]
    x_names <- gsub("d.", "", x_names)
    x_names <- substring(x_names, 1, unlist(gregexpr(".l", x_names))[1] - 1)
    if (s > 0) {
      x0 <- matrix(object[["x_x"]][, 1:m], NROW(object[["x_x"]]))
      x0 <- x0 + object[["w_x"]]
      x_diff <- matrix(object[["x_x"]][1, ], m)
      x_diff <- matrix(x_diff[, s:1], m)
      x_level <- x_diff * NA
      x_level[, s] <- x0[1, ] - x_diff[, s]
      if (s > 1) {
        for (i in (s - 1):1){
          x_level[, i] <- x_level[, i + 1] - x_diff[, i]
        } 
      }
      x <- rbind(t(x_level), x0)
    }
  }
  
  if (NROW(y) < NROW(x)) {
    y <- rbind(matrix(NA, NROW(x) - NROW(y), k), y)
  }
  if (NROW(y) > NROW(x)) {
    x <- rbind(matrix(NA, NROW(y) - NROW(x), m), x)
  }
  
  y <- stats::ts(y)
  x <- stats::ts(x)
  
  # Generate bvar data matrices
  temp <- y
  temp_names <- y_names
  if (!is.null(p)) {
    for (i in 1:p) {
      temp <- cbind(temp, stats::lag(y, -i))
      if (nchar(p) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_names <- c(temp_names, paste0(y_names, ".l", i_temp))
    } 
  }
  if (m > 0) {
    temp <- cbind(temp, x)
    if (nchar(s) > 2) {
      i_temp <- paste0(rep(0, nchar(s)), collapse = "")
    } else {
      i_temp <- paste0(rep(0, 2), collapse = "")
    }
    temp_names <- c(temp_names, paste0("exo.", x_names, ".l", i_temp))
    for (i in 1:s) {
      temp <- cbind(temp, stats::lag(x, -i))
      if (nchar(s) > 2) {
        i_temp <- paste0(c(rep(0, nchar(s) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_names <- c(temp_names, paste0("exo.", x_names, ".l", i_temp))
    }
  }
  
  temp <- stats::na.omit(temp)
  dimnames(temp)[[2]] <- temp_names
  y <- stats::ts(as.matrix(temp[, 1:k]), class = c("mts", "ts", "matrix"),
                 end = tsp_temp[2], frequency = tsp_temp[3])
  dimnames(y) <- list(NULL, y_names)
  
  x <- NULL
  if (p > 0 | m > 0) {
    x <- stats::ts(as.matrix(temp[, -(1:k)]), class = c("mts", "ts", "matrix"),
                   end = tsp_temp[2], frequency = tsp_temp[3])
    dimnames(x) <- list(NULL, temp_names[-(1:k)])
  }
  
  if (!is.null(object[["A0"]])) {
    if (tvp[["A0"]]) {
      n_a0 <- k * k
      A0 <- matrix(NA, n_a0 * tt, draws)
      for (i in 1:tt) {
        for (draw in 1:draws) {
          A0[(i - 1) * n_a0 + 1:n_a0, draw] <- object[["A0"]][[i]][draw, ]
        }
      }
    } else {
      A0 <- t(object[["A0"]]) 
    }
  } else {
    A0 <- NULL
  }
  
  n <- 0
  n_Pi_d <- 0
  C <- NULL
  x_det <- NULL
  x_det_names <- NULL
  if (!is.null(object[["Pi_d"]])) {
    if (tvp[["Pi_d"]]) {
      n <- NCOL(object[["Pi_d"]][[1]]) / k
    } else {
      n <- NCOL(object[["Pi_d"]]) / k
    }
    n_Pi_d <- n
    if (tvp[["Pi_d"]]) {
      n_pi <- k * n_Pi_d
      C <- matrix(NA, n_pi * tt, draws)
      for (i in 1:tt) {
        for (draw in 1:draws) {
          C[(i - 1) * n_pi + 1:n_pi, draw] <- object[["Pi_d"]][[i]][draw, ]
        }
      }
    } else {
      C <- t(object[["Pi_d"]]) 
    }
    if (!is.null(object[["w_d"]])) {
      x_det_names <- dimnames(object[["w_d"]])[[2]]
      x_det <- object[["w_d"]]
    }
  }
  
  n_c <- 0
  C_ur <- NULL
  if (!is.null(object[["C"]])) {
    if (tvp[["C"]]) {
      n_c <- NCOL(object[["C"]][[1]]) / k
    } else {
      n_c <- NCOL(object[["C"]]) / k 
    }
    n <- n + n_c
    
    if (tvp[["C"]]) {
      n_c_ur <- k * n_c
      C_ur <- matrix(NA, n_c_ur * tt, draws)
      for (i in 1:tt) {
        for (draw in 1:draws) {
          C_ur[(i - 1) * n_c_ur + 1:n_c_ur, draw] <- object[["C"]][[i]][draw, ]
        }
      }
    } else {
      C_ur <- t(object[["C"]]) 
    }
    
    if (!is.null(C)) {
      C <- rbind(C, C_ur)
    } else {
      C <- C_ur
    }
  }
  
  if (!is.null(object[["x_d"]]) & n_c > 0) {
    x_det_names <- c(x_det_names, dimnames(object[["x_d"]])[[2]]) 
    x_temp <- object[["x_d"]]
    if (is.null(x_det)) {
      x_det <- x_temp
    } else {
      x_det <- cbind(x_det, x_temp) 
    }
  }
  
  if (!is.null(x_det)) {
    x_names <- c(dimnames(x)[[2]], x_det_names)
    x <- cbind(x, x_det)
    dimnames(x)[[2]] <- x_names
  }
  
  if (!is.null(object[["Sigma"]])) {
    if (tvp[["Sigma"]]) {
      n_sigma <- k * k
      Sigma <- matrix(NA, n_sigma * tt, draws)
      for (i in 1:tt) {
        for (draw in 1:draws) {
          Sigma[(i - 1) * n_sigma + 1:n_sigma, draw] <- object[["Sigma"]][[i]][draw, ]
        }
      }
    } else {
      Sigma <- t(object[["Sigma"]]) 
    }
  } else {
    Sigma <- NULL
  }
  
  if (!is.null(object[["data"]])) {
    data <- object[["data"]]
  } else {
    data <- NULL
  }
  
  if (!is.null(object[["exogen"]])) {
    exogen <- object[["exogen"]]
  } else {
    exogen <- NULL
  }
  
  object <- bvar(data = data, exogen = exogen, y = y, x = x, A0 = A0, A = A, B = B, C = C, Sigma = Sigma)
  
  vars <- c("A", "B", "C", "Sigma", "A0")
  for (i in vars) {
    if(!is.null(object[[i]])) {
      if (object[["specifications"]][["tvp"]][[i]]) {
        for (j in 1:tt) {
          attr(object[[i]][[j]], "mcpar") <- specs
        }
      } else {
        attr(object[[i]], "mcpar") <- specs 
      }
    }
  }
  
  return(object)
}
