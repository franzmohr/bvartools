#' Reduced Rank Estimation of a Bayesian Vector Error Correction Model
#' 
#' Estimates a country-specific VECX model using the algorithm of Koop et al. (2010)
#' 
#' @param data a list containing the neccesary data to estimate a GVAR country model. Usually an object produced by \code{\link{country_models}}.
#' @param iterations an integer of MCMC draws including burn-in (defaults to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler (defaults to 5000).
#' These draws do not enter the computation of posterior moments, forecasts etc.
#' @param thin an integer specifying the thinning factor for the MCMC output. Defaults to 10, which means that the
#' forecast sequences contain only every tenth draw of the original sequence. Set \code{thin = 1} to obtain the full MCMC sequence.
#' 
#' @return Some object.
#' 
#' @export
bvecx_rr <- function(data, iterations = 50000, burnin = 5000, thin = 10){
  #### Generate model data and obtain model specifications ####
  estimation_data <- gen_vecx(data)
  y <- estimation_data$y
  ect <- estimation_data$ect
  x <- estimation_data$x
  t <- dim(y)[2]
  
  structural <- data$specs$structural
  tvp <- data$specs$tvp
  sv <- data$specs$sv
  if (sv) {
    structural <- TRUE
  }
  
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
  
  if (structural) {
    structural_est <- n > 1
    n_A0 <- n * (n - 1) / 2
  } else {
    structural_est <- FALSE
  }
  
  n_alpha <- n * r
  n_beta <- n_ect * r
  n_alpha_B <- n_alpha + n_B
  
  #### Set priors ####
  prior <- data$priors
  shrinkage <- prior$Shrinkage$type != "none"
  
  # Non-cointegration block
  B_mu_prior <- prior$B$constant[[1]]
  B_V_i_prior <- prior$B$constant[[2]]
  if (tvp[2]) {
    B_Q_df_post <- t + data$priors$B$tvp[[1]]
    B_Q_V_prior <- data$priors$B$tvp[[2]]
  }
  #if (shrinkage){
  #  shrinkage_type <- prior$Shrinkage$type
  #  if (shrinkage_type == "BVS"){
  #    B_restricted_variables <- prior$Shrinkage$spec$B[[1]]
  #    B_lpr_include_prior <- prior$Shrinkage$spec$B[[2]]
  #    B_lpr_exclude_prior <- prior$Shrinkage$spec$B[[3]]
  #  }
  #}
  
  # Cointegration block
  if (ecm) {
    if (tvp[1]) {
      # alpha_mu_prior <- prior$Pi$constant$mu
      # alpha_V_i_prior <- diag(1 / 10, n_alpha)
      # beta_mu_prior <- matrix(0, n_beta)
      # beta_V_i_prior <- diag(1 / 10, n_beta)
      # alpha_Q_df_post <- t + prior$Pi$tvp[[1]]
      # alpha_Q_V_prior <- prior$Pi$tvp[[2]]
      # rho <- prior$Pi$tvp[[3]][1]
    } else {
      if (r == 0) {
        v_i <- matrix(0, 1, 1)
        P_i <- matrix(0, 1, 1)
        G_i <- matrix(0, 1, 1)
      } else {
        v_i <- prior$Pi$constant$V_i$v_i
        P_i <- prior$Pi$constant$V_i$P_i
        G_i <- prior$Pi$constant$V_i$G_i
        strachan_prior <- G_i == "Omega_i"
        if (strachan_prior) {
          G_i_temp <- G_i
          rm(G_i)
        }
      }
    }
    alpha_B_mu_prior <- rbind(matrix(0, n_alpha), B_mu_prior)
    alpha_B_V_i_prior <- rbind(cbind(diag(0, n_alpha), matrix(0, n_alpha, n_B)), cbind(matrix(0, n_B, n_alpha), B_V_i_prior))
  }
  
  ## Covariance
   if (!structural) {
     Omega_df_post <- t + prior$Omega$constant[[1]] + r
     Omega_V_prior <- prior$Omega$constant[[2]]
   } else {
     A0_mu_prior <- prior$A0$constant[[1]]
     A0_V_i_prior <- prior$A0$constant[[2]]
  #   if (tvp[3]) {
  #     A0_Q_df_post <- t + prior$A0$tvp[[1]]
  #     A0_Q_V_prior <- prior$A0$tvp[[2]]
  #   }
  #   if (shrinkage){
  #     if (shrinkage_type == "BVS"){
  #       A0_lpr_include_prior <- prior$Shrinkage$spec$A0[[2]]
  #       A0_lpr_exclude_prior <- prior$Shrinkage$spec$A0[[3]]
  #     }
  #   }
  #   
     if (sv) {
       stoch_vol <- prior$Sigma$sv
     } else {
       Sigma_df_post <- t + prior$Omega$constant[[1]] + r
       Sigma_V_prior <- prior$Omega$constant[[2]] 
     } 
  }
  
  #### Initialise ####
  # Non-Cointegration block
  if (tvp[2]) {
    y_0 <- y * NA
    B_filter <- matrix(0, n_B, t + 1)
    B <- B_filter
    B_T <- diag(1, n_B)
    B_Q <- diag(1, n_B)
    B_res <- cbind(diag(B_Q_V_prior), B_Q_df_post - t)
    diag(B_Q) <- 1 / apply(B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    rm(B_res)
    B_a0 <- matrix(0, n_B)
    if (r == 0) {
      Z_B <- kronecker(t(x), diag(1, n))
    }
  }
  #if (shrinkage) {
  #  B_Gamma <- diag(1, B_n)
  #  Z_x <- kronecker(t(x), diag(1, n))
  #  Z_B <- Z_x
  #} else {
  #  Z_B <- kronecker(t(x), diag(1, n))
  #}
  
  # Cointegration block
  if (ecm) {
    ect_x <- rbind(matrix(NA, r, t), x)
    alpha_B_filter <- matrix(0, n_alpha_B, t + 1)
    if (tvp[1]) {
      #  rho <- .9999
      #  y_alpha_0 <- y * NA
      #  y_beta_0 <- y * NA
      #  alpha_filter <- matrix(0, n_alpha, t + 1)
      #  beta_filter <- matrix(0, n_beta, t + 1)
      #  Z_alpha <- matrix(NA, n * t, n_alpha)
      #  Z_beta <- matrix(NA, n * t, n_beta)
      #  #alpha_V_i_prior <- diag(1 / (1 - rho^2), n_alpha)
      #beta_V_prior <- diag(1 / (1 - rho^2), n_beta)
      #  alpha_Q <- diag(.0001, n_alpha)
      #  beta_Q <- diag(.0001, n_beta)
      #  Pi <- matrix(0, n_ect * n, t + 1)
      #  beta <- matrix(stats::rnorm(n_beta), n_ect * r, t + 1)
    } else {
      Pi <- matrix(0, n, n_ect)
      beta <- matrix(stats::rnorm(n_beta), n_ect, r)
    }
  }
  
  # Error block
  #y_dot <- y * 0
  
  ## Variance block
  if (structural) {
    #   y_dot_dot <- y * 0
    Sigma2 <- diag(apply(y, 1, stats::var), n)
    Sigma2_i <- solve(Sigma2)
    #   if (tvp[3]) {
    #     A0 <- t(matrix(diag(1, n), n, n * (t + 1)))
    #     A0_filter <- matrix(0, n_A0, t + 1)
    #     A0_Q <- diag(.0001, n_A0)
    #   } else {
    #     A0 <- diag(1, n)
    #   }
    #   if (shrinkage) {
    #     if (structural_est) {
    #       Gamma_A0 <- diag(1, n_A0)
    #     }
    #   }
    # }
    # 
    # if ((structural & tvp[3]) | sv) {
    #   Sigma2 <- t(matrix(diag(apply(y, 1, stats::var), n), n, n * t))
    #   Sigma2_i <- t(matrix(solve(diag(apply(y, 1, stats::var), n)), n, n * t))
    #   Omega <- Sigma2
    #   Omega_i <- Sigma2_i
  } else {
    Omega <- diag(apply(y, 1, stats::var), n)
    Omega_i <- solve(Omega)
  }
  
  #### Data containers ####
  # Non-cointegration block
  if (tvp[2]) {
    draws_B <- array(NA, c(n_B, t, store))
  } else {
    draws_B <- array(NA, c(n_B, 1, store))
  }
  
  if (shrinkage) {
    draws_B_shrinkage <- matrix(NA, n_B, iterations - burnin)
    if (structural_est) {
      draws_A0_shrinkage <- matrix(NA, n_A0, store)
    }
  }
  
  # Cointegration block
  if (tvp[1]) {
    draws_Pi <- array(NA, c(n * n_ect, t, store))
  } else {
    draws_Pi <- array(NA, c(n * n_ect, 1, store))
  }
  if (r == 0) {draws_Pi[,,] <- 0}
  
  # Covariance matrix
  if (sv | (tvp[3] & structural)) {
    draws_Omega <- array(NA, c(n * n, t, store))
  } else {
    draws_Omega <- array(NA, c(n * n, 1, store))
  }
  
  # Likelihood
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
    if (ecm) {
      if (tvp[2]) {
        if (tvp[1]) {
          for (i in 1:t) {
            ect_x[1:r, i] <- crossprod(matrix(beta[,i], n_ect), ect[, i]) 
          }
        } else {
          ect_x[1:r, ] <- crossprod(beta, ect) 
        }
        
        for (i in 1:t) {
          y_0[, i] <- y[, i] - matrix(alpha_B_filter[, i], n) %*% ect_x[, i]
        }
        Z_alpha_B <- kronecker(t(ect_x), diag(1, n))
        
        if (tvp[1]) {
          stop("Not implemented yet.")
        } else {
          if (strachan_prior) {
            G_i_temp <- Omega_i
          }
          alpha_B_V_i_prior[1:n_alpha, 1:n_alpha] <- kronecker(v_i * (crossprod(beta, P_i) %*% beta), G_i_temp) 
        }
        alpha_B_0 <- posterior_normal_sur(y_0, Z_alpha_B, Omega_i, alpha_B_mu_prior, alpha_B_V_i_prior)
        y_filter <- y - matrix(alpha_B_0, n) %*% ect_x
        
        if (tvp[1]) {
          #alpha_B_filter <- dk(y_filter, Z_alpha_B, Omega, B_Q, B_T, B_a0, B_Q) 
        } else {
          alpha_B_filter[-(1:n_alpha), ] <- dk(y_filter, Z_alpha_B[,-(1:n_alpha)], Omega, B_Q, B_T, B_a0, B_Q) 
        }
        B <- alpha_B_filter[-(1:n_alpha), ] + matrix(alpha_B_0[-(1:n_alpha), ], n_B, t + 1)
        
        B_res <- cbind(matrix(diag(B_Q_V_prior) + rowSums((B[,-1] - B[, -(t + 1)])^2)), B_Q_df_post)
        diag(B_Q) <- 1 / apply(B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
        
        y_hat <- y * NA
        for (i in 1:t) {
          y_hat[, i] <- y[, i]# - matrix(B[,i], n) %*%  x[, i]
        }
        
        if (tvp[1]) {
          stop("Not implemented yet.")
        } else {
          alpha <- matrix(alpha_B_0[1:n_alpha,], n)
          beta_temp <- posterior_koop2010_beta_sur(y_hat, ect, alpha, Omega_i, Omega_i, v_i, P_i)
          Pi <- beta_temp[[1]]
          beta <- beta_temp[[3]]
          y_dot <- beta_temp[[4]]
        }
      } else {
        if (sv | tvp[3]) {
          B_temp <- posterior_koop2010_sur(y, x, ect, r, beta, Omega_i, Omega_i, v_i, P_i, B_mu_prior, B_V_i_prior)
        } else {
          B_temp <- posterior_koop2010(y, x, ect, r, beta, Omega_i, Omega_i, v_i, P_i, B_mu_prior, B_V_i_prior)
        }
        B <- B_temp[[1]]
        if (r > 0) {
          Pi <- B_temp[[2]]
          beta <- B_temp[[4]]
        }
        y_dot <- B_temp[[5]]
      }
    } else {
      if (tvp[2]) {
        for (i in 1:t) {
          y_0[, i] <- y[, i] - matrix(B_filter[, i], n) %*% x[, i]
        }
        B_0 <- posterior_normal_sur(y_0, Z_B, Omega_i, B_mu_prior, B_V_i_prior)
        
        y_filter <- y - matrix(Z_B %*% B_0, n)
        B_filter <- dk(y_filter, Z_B, Omega, B_Q, B_T, B_a0, B_Q)
        B <- B_filter + matrix(B_0, n_B, t + 1)
      } else {
        if (sv | tvp[3]) {
          B <- matrix(posterior_normal_sur(y, Z_B, Omega_i, B_mu_prior, B_V_i_prior), n)
        } else {
          B <- matrix(posterior_normal(y, x, Omega_i, B_mu_prior, B_V_i_prior), n)
        }
        y_dot <- y - B %*% x
      }
    }
    
    if (structural) {
      if (structural_est) {
        for (i in 2:n) {
          p1 <- ((i - 1) * (i - 2) / 2 + 1) : ((i - 1) * i / 2)
          y_dot_temp <- matrix(y_dot[i,], 1)
          x_dot_temp <- matrix(-y_dot[1:(i - 1),], i - 1)
          if (sv) {
            Sigma2_i_temp <- matrix(Sigma2_i[(0:(t - 1) * n) + i, i])
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
            #if (shrinkage) {
            #  Gamma_A0[p1, p1] <- bvs(y_dot_temp, t(x_dot_temp), matrix(A0_temp[, -(t + 1)], i - 1), matrix(Gamma_A0[p1, p1], i - 1),
            #                          Sigma2_i_temp, matrix(1:(i - 1)), A0_lpr_include_prior, A0_lpr_exclude_prior)
            #  A0_temp <- Gamma_A0[p1, p1] %*% A0_temp
            #}
            
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
      }
      
      if (sv) {
        stop("Not implemented yet.")
      } else {
        stop("Not implemented yet.") 
      }
    } else {
      Omega_i <- rWishart(1, Omega_df_post, solve(tcrossprod(y_dot)))[,, 1]
      Omega <- solve(Omega_i)
    }
    
    #### Chunks ####
    #   # Korobilis (2013)
    #   if (shrinkage) {
    #     B_Gamma <- bvs(y_tilde, Z_x, B[, -(t + 1)], B_Gamma, Omega_i, B_restricted_variables, B_lpr_include_prior, B_lpr_exclude_prior)
    #     B <- B_Gamma %*% B
    # #### Error covariance matrix ####
    # if (structural) {
    #   #### Structural ####
    
    #     
    #     #### Volatility state ####
    #     if (tvp[3]) {
    #       for (i in 1:t) {
    #         y_dot_dot[, i] <- A0[(i - 1) * n + 1:n, ] %*% y_dot[, i]
    #       }
    #     } else {
    #       y_dot_dot <- A0 %*% y_dot
    #     }
    #     
    #     if (sv) {
    #       for (i in 1:n) {
    #         sv_temp <- stochvol::svsample2(y_dot_dot[i,] + .0000001,
    #                                        priormu = stoch_vol[[i]]$priors$priormu,
    #                                        priorphi = stoch_vol[[i]]$priors$priorphi,
    #                                        priorsigma = stoch_vol[[i]]$priors$priorsigma,
    #                                        startpara = stoch_vol[[i]]$para,
    #                                        startlatent = stoch_vol[[i]]$latent)
    #         stoch_vol[[i]]$para <- sv_temp$para
    #         stoch_vol[[i]]$latent <- sv_temp$latent
    #         h <- exp(sv_temp$latent)
    #         Sigma2[i + 0:(t - 1) * n, i] <- h
    #         Sigma2_i[i + 0:(t - 1) * n, i] <- 1 / h
    #       }
    #     } else {
    #       Sigma_res <- cbind(diag(Sigma_V_prior) + rowSums(y_dot_dot^2), Sigma_df_post)
    #       diag(Sigma2_i) <- apply(Sigma_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    #       Sigma2 <- solve(Sigma2_i)
    #     }
    #     
    #     if (tvp[3]) {
    #       for (i in 1:t) {
    #         A0_i_temp <- solve(A0[(i - 1) * n + 1:n, ])
    #         if (sv) {
    #           Omega[(i - 1) * n + 1:n,] <- A0_i_temp %*% tcrossprod(Sigma2[(i - 1) * n + 1:n,], A0_i_temp) 
    #         } else {
    #           Omega[(i - 1) * n + 1:n,] <- A0_i %*% tcrossprod(Sigma2, A0_i)
    #         }
    #         Omega_i[(i - 1) * n + 1:n,] <- solve(Omega[(i - 1) * n + 1:n,])
    #       }
    #     } else {
    #       A0_i <- solve(A0)
    #       if (sv) {
    #         for (i in 1:t) {
    #           Omega[(i - 1) * n + 1:n, ] <- A0_i %*% tcrossprod(Sigma2[(i - 1) * n + 1:n, ], A0_i)
    #           Omega_i[(i - 1) * n + 1:n, ] <- solve(Omega[(i - 1) * n + 1:n, ]) 
    #         }
    #       } else {
    #         Omega <- A0_i %*% tcrossprod(Sigma2, A0_i)
    #         Omega_i <- solve(Omega) 
    #       }
    #     }
    #   } else {
    #     if (sv) {
    #       sv_temp <- stochvol::svsample2(y_dot + .0000001,
    #                                      priormu = stoch_vol[[1]]$priors$priormu,
    #                                      priorphi = stoch_vol[[1]]$priors$priorphi,
    #                                      priorsigma = stoch_vol[[1]]$priors$priorsigma,
    #                                      startpara = stoch_vol[[1]]$para,
    #                                      startlatent = stoch_vol[[1]]$latent)
    #       stoch_vol[[1]]$para <- sv_temp$para
    #       stoch_vol[[1]]$latent <- sv_temp$latent
    #       h <- exp(sv_temp$latent)
    #       Sigma2[i + 0:(t - 1) * n, i] <- h
    #       Sigma2_i[i + 0:(t - 1) * n, i] <- 1 / h
    #     } else {
    #       Omega_i <- matrix(stats::rgamma(1, shape = Sigma_df_post / 2, rate = (Sigma_V_prior  + sum(y_dot^2)) / 2))
    #       Omega <- 1 / Omega_i 
    #     }
    #   }
    
    # 
    # #### Post-burnin procedures ####
    if (draw > burnin){
      pos_draw <- draw - burnin
      if (tvp[2]) {
        draws_B[, , pos_draw] <- B[, -(t + 1)]
      } else {
        draws_B[, 1, pos_draw] <- B
      }
      if (shrinkage) {
        draws_B_shrinkage[, pos_draw] <- diag(B_Gamma)
        if (structural_est) {
          draws_A0_shrinkage[, pos_draw] <- diag(Gamma_A0)
        }
      }
      if (ecm){
        if (tvp[1]) {
          draws_Pi[,, pos_draw] <- Pi[, -(t + 1)] 
        } else {
          draws_Pi[, 1, pos_draw] <- Pi
        }
      }
      
      if (sv | (tvp[3] & structural)) {
        for (i in 1:t) {
          draws_Omega[, i, pos_draw] <- Omega[(i - 1) * n + 1:n,] 
        }
      } else {
        draws_Omega[, 1, pos_draw] <- Omega
      }
      
      draws_LL[, pos_draw] <- getLL(y_dot, Omega, Omega_i)
      ###########################################################################
      if (FALSE) {
        p <- estimation_data$domestic["lag"]
        graphics::par(mfcol = c(3, 3))
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter y")
        #stats::plot.ts(cbind(0, draws_B[n * n * p + 1, 1, ]), plot.type = "single", col = "black", ylab = "Contemp y")
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic y")
        #stats::plot.ts(cbind(0, draws_B[n_B - 2, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic y")
        
        stats::plot.ts(cbind(0, draws_Pi[n_ect - 2, 1, ]), plot.type = "single", col = "black", ylab = "Pi y")
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + n + 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter Dp")
        #stats::plot.ts(cbind(0, draws_B[n * n * p + n + 2, 1, ]), plot.type = "single", col = "black", ylab = "Contemp Dp")
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic Dp")
        #stats::plot.ts(cbind(0, draws_B[n_B - 1, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic Dp")
        
        stats::plot.ts(cbind(0, draws_Pi[n_ect - 1, 1, ]), plot.type = "single", col = "black", ylab = "Pi Dp")
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 2 * n + 3, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter r")
        #stats::plot.ts(cbind(0, draws_B[n * n * p + 2 * n + 3, 1, ]), plot.type = "single", col = "black", ylab = "Contemp r")
        
        stats::plot.ts(cbind(0, t(apply(draws_B[n_B, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic r")
        #stats::plot.ts(cbind(0, draws_B[n_B, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic r")
        
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
  
  B_t <- ifelse(tvp[2], t, 1)
  if (tvp[2]) {
    B <- draws_B[,, tw]
  } else {
    B <- array(draws_B[,,tw], dim = c(n_B, 1, draws))
  }
  
  Pi_t <- ifelse(tvp[1], t, 1)
  if (tvp[1]) {
    Pi <- draws_Pi[,, tw]
  } else {
    Pi <- array(draws_Pi[,, tw], dim = c(n * n_ect, 1, draws)) 
  }
  
  if (sv | (structural & tvp[3])) {
    Omega <- draws_Omega[,, tw]
  } else {
    Omega <- array(draws_Omega[, 1, tw], dim = c(n * n, 1, draws))
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