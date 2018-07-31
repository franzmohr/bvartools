#' Estimate a Bayesian Vector Error Correction Model
#' 
#' Estimates a country-specific VECX model using the algorithm of Koop et al. (2010)
#' 
#' @param data a list containing the neccesary data to estimate a GVAR country model.
#' @param iterations an integer of MCMC draws including burn-in (defaults to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler (defaults to 5000).
#' These draws do not enter the computation of posterior moments, forecasts etc.
#' @param thin an integer specifying the thinning factor for the MCMC output. Defaults to 10, which means that the
#' forecast sequences contain only every tenth draw of the original sequence. Set \code{thin = 1} to obtain the full MCMC sequence.
#' 
#' @return Some object.
#' 
#' @export
vecx2varx <- function(data){
  #### Transform VECX parameters to VARX ####
  names_y <- gsub("d.", "", dimnames(y)[[1]])
  names_ect <- dimnames(ect)[[1]]
  names_x <- gsub("d.", "", dimnames(x)[[1]])
  
  if (is.null(thin)){
    tw <- seq(to = store, length.out = store)
  } else {
    tw <- sample(x = 1:store, size = store / thin, replace = FALSE)
  }
  
  draws <- length(tw)
  A0_t <- ifelse(tvp[3], t, 1)
  A0 <- array(diag(1, n), dim = c(n * n, A0_t, draws))
  dimnames(A0) <- list(paste(names_y, rep(names_y, each = n), sep = "_"), NULL, NULL)
  
  B_t <- ifelse(tvp[1], t, 1)
  if (tvp[1]) {
    B <- draws_B[,, tw]
  } else {
    B <- array(draws_B[, tw], dim = c(n_B, 1, draws))
  }
  
  Pi_t <- ifelse(tvp[2], t, 1)
  if (tvp[2]) {
    Pi <- draws_Pi[,, tw]
  } else {
    Pi <- array(draws_Pi[, tw], dim = c(n * n_ect, 1, draws)) 
  }
  
  if (sv | (structural & tvp[3])) {
    Omega <- draws_Omega[,, tw]
  } else {
    Omega <- array(draws_Omega[, tw], dim = c(n * n, 1, draws))
  }
  dimnames(Omega) <- list(paste(names_y, rep(names_y, each = n), sep = "_"), NULL, NULL)
  
  # Domestic variables
  p <- estimation_data$domestic["lag"]
  n_d <- estimation_data$domestic["total"]
  if (p > 0){
    B_d <- array(B[1:(n * n_d),,], c(n * n_d, B_t, draws))
  }
  
  A_d_t <- max(B_t, Pi_t)
  p <- p + 1
  if (p > 1){
    W <- diag(-1, n * p)
    W[1:n, 1:n] <- diag(1, n)
    W[-(1:n), -(n * (p - 1) + 1:n)] <- W[-(1:n), -(n * (p - 1) + 1:n)] + diag(n * (p - 1))
    J <- matrix(0, n, n * p)
    J[1:n, 1:n] <- diag(1, n)
    
    A_d <- array(NA, c(n * n * p, A_d_t, draws))
    for (draw in 1:draws){
      for (i in 1:A_d_t) {
        A_d[, i, draw] <- cbind(matrix(Pi[1:(n * n), ifelse(Pi_t == 1, 1, i), draw], n),
                                matrix(B_d[, ifelse(B_t == 1, 1, i), draw], n)) %*% W + J 
      }
    }
    names_d <- NULL
    for (i in 1:p){
      names_d <- c(names_d, paste(names_y, ".l", i, sep = ""))
    }
    dimnames(A_d) <- list(paste(names_y, rep(names_d, each = n), sep = "_"), NULL, NULL)
  } else {
    A_d <- array(NA, c(n * n, Pi_t, draws))
    for (i in 1:Pi_t) {
      A_d[,i,] <- Pi[1:(n * n), i, ] + matrix(diag(1, n), n * n, draws)
    }
    dimnames(A_d) <- list(paste(names_y, rep(names_y, each = n), sep = "_"), NULL, NULL)
  }
  
  # Foreign variables
  n_star <- estimation_data$foreign["dim"]
  p_star <- estimation_data$foreign["lag"]
  n_s <- estimation_data$foreign["total"]
  B_s <- array(B[n_d * n + 1:(n * n_s),,], c(n * n_s, B_t, draws))
  
  p_star <- p_star + 1
  if (p_star > 0){
    W.s <- diag(-1, n_star * (p_star + 1))
    W.s[1:n_star, 1:n_star] <- 0
    W.s[1:n_star, n_star + 1:n_star] <- diag(1, n_star)
    W.s[-(1:n_star), 1:(n_star * p_star)] <- W.s[-(1:n_star), 1:(n_star * p_star)] + diag(1, n_star * p_star)
    
    A_s <- array(NA, c(n * n_star * (p_star + 1), A_d_t, draws))
    for (draw in 1:draws){
      for (i in 1:A_d_t) {
        A_s[, i, draw] <- cbind(matrix(Pi[n * n + 1:(n * n_star), ifelse(Pi_t == 1, 1, i), draw], n),
                                matrix(B_s[, ifelse(B_t == 1, 1, i), draw], n)) %*% W.s
      }
    }
    
    names_s <- names_x[n_d + 1:n_star]
    for (i in 1:p_star){
      names_s <- c(names_s, paste(names_x[n_d + 1:n_star], ".l", i, sep = ""))
    }
    dimnames(A_s) <- list(paste(names_y, rep(names_s, each = n), sep = "_"), NULL, NULL)
  } else {
    stop("Lag of foreign star variables in VARX model must be at least 1.")
  }
  
  A_s_0 <- array(A_s[1:(n_star * n), ,], c(n * n_star, B_t, draws))
  dimnames(A_s_0) <- list(dimnames(A_s)[[1]][1:(n * n_star)], NULL, NULL)
  names_A_s <- dimnames(A_s)[[1]][-(1:(n * n_star))]
  A_s <- array(A_s[-(1:(n * n_star)), ,], c(n * n_s, B_t, draws))
  dimnames(A_s) <- list(names_A_s, NULL, NULL)
  
  # Global variables
  A_g <- NULL
  global <- !is.null(data$x.g)
  n_global <- estimation_data$global["dim"]
  p_global <- estimation_data$global["lag"]
  n_g <- estimation_data$global["total"]
  if (global){
    if (!is.na(p_global)){
      names_g <- gsub("d.", "", dimnames(x)[[1]][n_d + n_s + 1:n_global])
      
      B_g <- array(B[n * n_d + n * n_s + 1:(n * n_g), ,], c(n * n_g, B_t, draws))
      A_g <- array(NA, c(n * (n_g + 1), A_d_t, draws))
      if (p_global > 1){
        W_g <- diag(-1, n_global * (p_global + 1))
        W_g[1:n_global, 1:n_global] <- 0
        W_g[1:n_global, n_global + 1:n_global] <- diag(1, n_global)
        W_g[-(1:n_global), 1:(n_global * p_global)] <- W_g[-(1:n_global), 1:(n_global * p_global)] + diag(1, n_global * p_global)
        
        for (draw in 1:draws){
          for (i in 1:A_d_t) {
            A_g[, i, draw] <- cbind(matrix(Pi[n * (n + n_star) + 1:(n * n_global), ifelse(Pi_t == 1, 1, i), draw], n), 
                                    matrix(B_g[, ifelse(B_t == 1, 1, i), draw], n)) %*% W_g
          }
        }
        
        names_global <- names_g
        for (i in 1:p_global) {
          names_global <- c(names_global, paste(names_g, ".l", i, sep = ""))
        }
        dimnames(A_g) <- list(paste(names_y, rep(names_global, each = n), sep="_"), NULL, NULL)
      } else {
        
        for (draw in 1:draws){
          for (i in 1:A_d_t) {
            A_g[, i, draw] <- cbind(matrix(B_g[, ifelse(B_t == 1, 1, i), draw], n),
                                    matrix(Pi[n * (n + n_star) + 1:(n * n_global), ifelse(Pi_t == 1, 1, i), draw], n) - matrix(B_g[, ifelse(B_t == 1, 1, i), draw], n)) 
          }
        }
        names_global <- c(names_x[n_d + n_s + 1:n_global], paste(names_x[n_d + n_s + 1:n_global], ".l", 1, sep = ""))
        dimnames(A_g) <- list(paste(names_y, rep(names_global, each = n), sep = "_"), NULL, NULL)
      }
    }
  }
  
  # Deterministic terms
  n_det <- n_det_r + n_det_ur
  if (n_det > 0) {
    det <- array(NA, c(n_det * n, A_d_t, draws))
    names_det <- NULL
    # Unrestricted determinisitc terms
    if (n_det_ur > 0){
      for (i in 1:A_d_t) {
        det[1:(n * n_det_ur), i,] <- B[n * (n_d + n_s + n_g) + 1:(n_det_ur * n), ifelse(B_t == 1, 1, i),]
      }
      names_det <- c(names_det, names_x[n_d + n_s + n_g + 1:n_det_ur]) 
    }
    # Restriced deterministic terms
    if (n_det_r > 0){
      for (i in 1:A_d_t) {
        det[(n * n_det_ur) + 1:(n * n_det_r), i,] <- Pi[n * (n + n_star + n_global) + 1:(n * n_det_r), ifelse(Pi_t == 1, 1, i),] 
      }
      names_det <- c(names_det, names_ect[n + n_star + n_global + 1:n_det_r])
    }
    dimnames(det) <- list(paste(names_y, rep(names_det, each = n), sep = "_"), NULL, NULL)
  } else {
    det <- NULL
  }
  
  #### Test statistics ####
  draws_LL <- draws_LL[, tw]
  test_stat <- c("AIC" = 2 * sum(log(rowMeans(draws_LL))) - (n_alpha_B) * 2,
                 "BIC" = 2 * sum(log(rowMeans(draws_LL))) - (n_alpha_B) * log(t),
                 "HQ" = 2 * sum(log(rowMeans(draws_LL))) - (n_alpha_B) * 2 * log(log(t)))
  
  ##### Prepare output object ####
  data$specs$tvp <- (tvp[1] | tvp[2])
  result <- c(data, list(coefs = list(A0 = A0, A_d = A_d, A_s_0 = A_s_0, A_s = A_s, A_g = A_g, A_det = det, Omega = Omega),
                         criteria = list("Criteria" = test_stat)))
  
  if (shrinkage) {
    A_p_include <- matrix(rowMeans(draws_B_shrinkage), n)
    dimnames(A_p_include) <- list(names_y, names_x)
    result$coefs$A_p_include <- A_p_include
    
    if (structural_est) {
      A0_p_include <- matrix(NA, n, n)
      A0_means <- rowMeans(draws_A0_shrinkage)
      for (i in 2:n) {
        A0_p_include[i, 1:(i-1)] <- A0_means[((i - 1) * (i - 2) / 2 + 1) : ((i - 1) * i / 2)]
      }
      dimnames(A0_p_include) <- list(names_y, names_y)
      result$coefs$A0_p_include <- A0_p_include
    }
  }
  return(result)
}