#' Estimate a Vector Error Correction Model
#' 
#' Estimates a country-specific VECX model.
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
bvecx_2 <- function(data, iterations = 15000, burnin = 5000, thin = NULL){
  #### Generate model data and obtain model specifications ####
  estimation_data <- gen_vecx(data)
  y <- estimation_data$y
  ect <- estimation_data$ect
  x <- estimation_data$x
  t <- dim(y)[2]
  
  r <- estimation_data$r
  if (is.na(r)) {stop("Rank of Pi must be specified.")}
  if (r != 0) {ecm <- TRUE} else {ecm <- FALSE}
  
  store <- iterations - burnin
  if (store > (iterations - burnin)) {stop("Number of iterations must be larger than burn-in draws.")}
  
  n <- dim(y)[1]
  n_ect <- dim(ect)[1]
  n_x <- dim(x)[1]
  n_B <- n * n_x
  n_det_r <- estimation_data$deterministic["restricted"]
  n_det_ur <- estimation_data$deterministic["unrestricted"]
  
  structural_est <- n > 1
  n_A0 <- n * (n - 1) / 2
  
  n_alpha <- n * r
  n_beta <- n_ect * r
  
  n_alpha_B <- n_alpha + n_B
  
  #### Set priors ####
  prior <- data$priors
  shrinkage <- prior$Shrinkage$type != "none"
  
  # Non-cointegration block
  B_mu_prior <- prior$B$constant[[1]]
  B_V_i_prior <- prior$B$constant[[2]]
  B_Q_df_post <- t + data$priors$B$tvp[[1]]
  B_Q_V_prior <- data$priors$B$tvp[[2]]
  
  # Cointegration block
  v_i <- prior$Pi$constant$V_i$v_i
  P_i <- prior$Pi$constant$V_i$P_i

  # Covariance matrix
  ## Constant covariance
  if (FALSE) {
  if (!structural) {
    Omega_df_post <- t + prior$Omega$constant[[1]]
    if (r != "full") {
      Omega_df_post <- Omega_df_post + r
    }
    Omega_V_prior <- prior$Omega$constant[[2]]
  } else {
    A0_mu_prior <- prior$A0$constant[[1]]
    A0_V_i_prior <- prior$A0$constant[[2]]
    if (tvp[3]) {
      A0_Q_df_post <- t + prior$A0$tvp[[1]]
      A0_Q_V_prior <- prior$A0$tvp[[2]]
    }
    if (shrinkage){
      if (shrinkage_type == "BVS"){
        A0_lpr_include_prior <- prior$Shrinkage$spec$A0[[2]]
        A0_lpr_exclude_prior <- prior$Shrinkage$spec$A0[[3]]
      }
    }
    
    if (sv) {
      stoch_vol <- prior$Sigma$sv
    } else {
      Sigma_df_post <- t + prior$Omega$constant[[1]]
      Sigma_V_prior <- prior$Omega$constant[[2]] 
    } 
  }
  }
  
  #### Initialise ####
  # Non-Cointegration block
  y_tilde <- y * NA
  y_tilde_0 <- y * NA
  ect_x <- rbind(ect, x)
  #y_tilde_filter <- y * NA
  Z_B <- kronecker(t(x), diag(1, n))
  #ect_x <- rbind(matrix(NA, r, t), x)
  B_filter <- matrix(0, n_B, t + 1)
  B <- B_filter
  B_T <- diag(1, n_B)
  B_Q <- diag(1, n_B)
  B_a0 <- matrix(0, n_B)
  B_res <- cbind(diag(B_Q_V_prior), B_Q_df_post - t)
  diag(B_Q) <- 1 / apply(B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
  rm(B_res)
  
  # Cointegration block
  y_hat <- y * NA 
  #y_beta_0 <- y * NA
  #beta_filter <- matrix(0, n_beta, t + 1)
  #Z_beta <- matrix(NA, n * t, n_beta)
  #beta_Q <- diag(1, n_beta)
  Pi <- matrix(0, n, n_ect)
  beta <- matrix(stats::rnorm(n_beta), n_ect, r) 

  # Error block
  #y_dot <- y * 0
  Omega <- diag(apply(y, 1, stats::var), n)
  Omega_i <- solve(Omega)

  #### Data containers ####
  draws_B <- array(NA, c(n_B, t, store))
  draws_Pi <- array(NA, c(n * n_ect, 1, store))
  draws_Omega <- array(NA, c(n * n, 1, store))
  draws_LL <- matrix(NA, t, store)
  
  # Progress bar
  if (iterations <= 100){
    update <- 1:iterations
  } else {
    update <- round(seq(from = 0, to = iterations, length.out = 100)) # Refreshment rate of the progress bar
  }
  pb <- utils::txtProgressBar(width = 70, style = 3)
  #### Start sampling ####
  for (draw in 1:iterations){
    B_temp <- posterior_koop2010(y, x, ect, beta, Omega_i, Omega_i, v_i, P_i, B_mu_prior, B_V_i_prior)
    
    B <- B_temp[[1]]
    Pi <- B_temp[[2]]
    beta <- B_temp[[4]]
    
    y_dot <- y - Pi %*% ect - B %*% x
    
    
    Omega_i <- rWishart(1, t + r, solve(tcrossprod(y_dot)))[,, 1]
    Omega <- solve(Omega_i)
    
  if (FALSE)  { 
    #### Error covariance matrix ####
    if (structural) {
      #### Structural ####
      if (structural_est) {
        for (i in 2:n) {
          p1 <- ((i - 1) * (i - 2) / 2 + 1) : ((i - 1) * i / 2)
          y_dot_temp <- matrix(y_dot[i,], 1)
          x_dot_temp <- matrix(-y_dot[1:(i - 1),], i - 1)
          if (sv) {
            Sigma2_i_temp <- matrix(Sigma2_i[(0:(t - 1) * n) + i,i])
          } else {
            Sigma2_i_temp <- matrix(diag(Sigma2_i)[i]) 
          }
          A0_mu_prior_temp <- matrix(A0_mu_prior[p1,], i - 1)
          A0_V_i_prior_temp <- matrix(A0_V_i_prior[p1,p1], i - 1)
          
          if (tvp[3]) {
            y_dot_temp_0 <- y_dot_temp * 0
            A0_filter_temp <- matrix(A0_filter[p1,], i - 1)
            if (sv) {
              Sigma2_temp <- matrix(Sigma2[(0:(t - 1) * n) + i,i])
            } else {
              Sigma2_temp <- matrix(diag(Sigma2)[i]) 
            }
            
            # Initialise
            for (j in 1:t) {
              y_dot_temp_0[, j] <- y_dot_temp[, j] - x_dot_temp[, j] %*% A0_filter_temp[, j]
            }
            A0_0 <- matrix(posterior_normal_sur(y_dot_temp_0, t(x_dot_temp), Sigma2_i_temp, A0_mu_prior_temp, A0_V_i_prior_temp), 1)
            y_dot_temp_filter <- y_dot_temp - A0_0 %*% x_dot_temp
            
            # Durbin & Koopman (2002)
            A0_filter_temp <- dk(y_dot_temp_filter, t(x_dot_temp), Sigma2_temp, matrix(A0_Q[p1, p1], i - 1), diag(1, i - 1),
                                 matrix(0, i - 1), matrix(A0_Q[p1, p1], i - 1))
            A0_temp <- A0_filter_temp + matrix(A0_0, i - 1, t + 1)
            
            # Korobilis (2013)
            if (shrinkage) {
              Gamma_A0[p1, p1] <- bvs(y_dot_temp, t(x_dot_temp), matrix(A0_temp[, -(t + 1)], i - 1), matrix(Gamma_A0[p1, p1], i - 1),
                                      Sigma2_i_temp, matrix(1:(i - 1)), A0_lpr_include_prior, A0_lpr_exclude_prior)
              A0_temp <- Gamma_A0[p1, p1] %*% A0_temp
            }
            
            # Covariance matrix of the coefficients
            A0_res <- cbind(matrix(diag(A0_Q_V_prior)[p1] + rowSums(matrix(A0_temp[,-1] - A0_temp[, -(t + 1)], i - 1)^2)), A0_Q_df_post[p1])
            diag(A0_Q)[p1] <- 1 / apply(A0_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
            
            A0[i + (0:t * n), 1:(i - 1)] <- t(A0_temp)
          } else {
            if (sv) {
              A0_temp <- posterior_normal_sur(y_dot_temp, t(x_dot_temp), Sigma2_i_temp, A0_mu_prior_temp, A0_V_i_prior_temp)
            } else {
              A0_temp <- posterior_normal(y_dot_temp, x_dot_temp, Sigma2_i_temp, A0_mu_prior_temp, A0_V_i_prior_temp) 
            }
            A0[i, 1:(i - 1)] <- A0_temp 
          }
        }
        
        #### Volatility state ####
        if (tvp[3]) {
          for (i in 1:t) {
            y_dot_dot[, i] <- A0[(i - 1) * n + 1:n, ] %*% y_dot[, i]
          }
        } else {
          y_dot_dot <- A0 %*% y_dot
        }
        
        if (sv) {
          for (i in 1:n) {
            sv_temp <- stochvol::svsample2(y_dot_dot[i,] + .0000001,
                                           priormu = stoch_vol[[i]]$priors$priormu,
                                           priorphi = stoch_vol[[i]]$priors$priorphi,
                                           priorsigma = stoch_vol[[i]]$priors$priorsigma,
                                           startpara = stoch_vol[[i]]$para,
                                           startlatent = stoch_vol[[i]]$latent)
            stoch_vol[[i]]$para <- sv_temp$para
            stoch_vol[[i]]$latent <- sv_temp$latent
            h <- exp(sv_temp$latent)
            Sigma2[i + 0:(t - 1) * n, i] <- h
            Sigma2_i[i + 0:(t - 1) * n, i] <- 1 / h
          }
        } else {
          Sigma_res <- cbind(diag(Sigma_V_prior) + rowSums(y_dot_dot^2), Sigma_df_post)
          diag(Sigma2_i) <- apply(Sigma_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
          Sigma2 <- solve(Sigma2_i)
        }
        
        if (tvp[3]) {
          for (i in 1:t) {
            A0_i_temp <- solve(A0[(i - 1) * n + 1:n, ])
            if (sv) {
              Omega[(i - 1) * n + 1:n,] <- A0_i_temp %*% tcrossprod(Sigma2[(i - 1) * n + 1:n,], A0_i_temp) 
            } else {
              Omega[(i - 1) * n + 1:n,] <- A0_i %*% tcrossprod(Sigma2, A0_i)
            }
            Omega_i[(i - 1) * n + 1:n,] <- solve(Omega[(i - 1) * n + 1:n,])
          }
        } else {
          A0_i <- solve(A0)
          if (sv) {
            for (i in 1:t) {
              Omega[(i - 1) * n + 1:n, ] <- A0_i %*% tcrossprod(Sigma2[(i - 1) * n + 1:n, ], A0_i)
              Omega_i[(i - 1) * n + 1:n, ] <- solve(Omega[(i - 1) * n + 1:n, ]) 
            }
          } else {
            Omega <- A0_i %*% tcrossprod(Sigma2, A0_i)
            Omega_i <- solve(Omega) 
          }
        }
      } else {
        if (sv) {
          sv_temp <- stochvol::svsample2(y_dot + .0000001,
                                         priormu = stoch_vol[[1]]$priors$priormu,
                                         priorphi = stoch_vol[[1]]$priors$priorphi,
                                         priorsigma = stoch_vol[[1]]$priors$priorsigma,
                                         startpara = stoch_vol[[1]]$para,
                                         startlatent = stoch_vol[[1]]$latent)
          stoch_vol[[1]]$para <- sv_temp$para
          stoch_vol[[1]]$latent <- sv_temp$latent
          h <- exp(sv_temp$latent)
          Sigma2[i + 0:(t - 1) * n, i] <- h
          Sigma2_i[i + 0:(t - 1) * n, i] <- 1 / h
        } else {
          Omega_i <- matrix(stats::rgamma(1, shape = Sigma_df_post / 2, rate = (Sigma_V_prior  + sum(y_dot^2)) / 2))
          Omega <- 1 / Omega_i 
        }
      }
    } else {
      Omega_i <- stats::rWishart(1, Omega_df_post, solve(tcrossprod(y_dot)))[,,1]
      Omega <- solve(Omega_i)  
    }
    
    
  }
    #### Post-burnin procedures ####
    if (draw > burnin){
      pos_draw <- draw - burnin
      draws_B[, 1,pos_draw] <- B#[,-(t + 1)]
      draws_Pi[, 1, pos_draw] <- Pi
      draws_Omega[, 1,pos_draw] <- Omega

      draws_LL[, pos_draw] <- getLL(y_dot, Omega, Omega_i)
      ###########################################################################
      if (FALSE) {
        p <- estimation_data$domestic["lag"]
        #if (draw %in% update & pos_draw > 10) {
          graphics::par(mfcol = c(3, 3))
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter y")
          stats::plot.ts(cbind(0, draws_B[n * n * p + 1, 1, ]), plot.type = "single", col = "black", ylab = "Contemp y")
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic y")
          stats::plot.ts(cbind(0, draws_B[n_B - 2, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic y")
          
          stats::plot.ts(cbind(0, draws_Pi[n_ect - 2, 1, ]), plot.type = "single", col = "black", ylab = "Pi y")
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + n + 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter Dp")
          stats::plot.ts(cbind(0, draws_B[n * n * p + n + 2, 1, ]), plot.type = "single", col = "black", ylab = "Contemp Dp")
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic Dp")
          stats::plot.ts(cbind(0, draws_B[n_B - 1, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic Dp")
          
          stats::plot.ts(cbind(0, draws_Pi[n_ect - 1, 1, ]), plot.type = "single", col = "black", ylab = "Pi Dp")
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 2 * n + 3, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter r")
          stats::plot.ts(cbind(0, draws_B[n * n * p + 2 * n + 3, 1, ]), plot.type = "single", col = "black", ylab = "Contemp r")
          
          #stats::plot.ts(cbind(0, t(apply(draws_B[n_B, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic r")
          stats::plot.ts(cbind(0, draws_B[n_B, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic r")
          
          stats::plot.ts(cbind(0, draws_Pi[n_ect, 1, ]), plot.type = "single", col = "black", ylab = "Pi r")
          #stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter")
          
          #for (i in 1:12) {
          #  stats::plot.ts(cbind(0, t(apply(draws_Pi[i, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #                 plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Cointegration parameter")
          #}

          #if (sv | tvp[3]) {
          #  stats::plot.ts(log(t(apply(draws_Omega[1, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
          #          plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Log-Volatility") 
          #}
          graphics::par(mfcol = c(1, 1))
        #}
      }
      ###########################################################################
    }
    # Update progress bar
    if (draw %in% update) {utils::setTxtProgressBar(pb, draw / iterations)}
  }
  close(pb)
  
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
    B <- array(draws_B[, tw], dim = c(B_n, 1, draws))
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
