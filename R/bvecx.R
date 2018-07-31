#' Estimation of a Bayesian Vector Error Correction Model
#' 
#' Estimates a country-specific VECX model without the cointegration matrix \eqn{\beta}.
#' 
#' @param data a list containing the neccesary data to estimate a GVAR country model. Usually and object produced by \code{\link{country_models}}.
#' @param iterations an integer of MCMC draws including burn-in (defaults to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler (defaults to 5000).
#' These draws do not enter the computation of posterior moments, forecasts etc.
#' @param thin an integer specifying the thinning factor for the MCMC output. Set \code{thin = 1} to obtain the full MCMC sequence.
#' 
#' @return Some object.
#' 
#' @export
bvecx <- function(data, iterations = 50000, burnin = 5000, thin = 10){
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
  
  store <- iterations - burnin
  if (store > (iterations - burnin)) {stop("Number of iterations must be larger than burn-in draws.")}
  
  n <- dim(y)[1]
  n_ect <- dim(ect)[1]
  n_x <- dim(x)[1]
  n_B <- n * n_x
  n_Pi <- n * n_ect
  n_Pi_B <- n_Pi + n_B
  n_det_r <- estimation_data$deterministic["restricted"]
  n_det_ur <- estimation_data$deterministic["unrestricted"]
  
  if (structural) {
    structural_est <- n > 1
    n_A0 <- n * (n - 1) / 2
  } else {
    structural_est <- FALSE
  }
  tvp_covar <- sv | (tvp[3] & structural)

  #### Set priors and fixed posteriors ####
  prior <- data$priors
  shrinkage <- prior$Shrinkage$type != "none"
  
  # Non-cointegration block
  Pi_B_mu_prior <- rbind(prior$Pi$constant$mu, prior$B$constant[[1]])
  Pi_B_V_i_prior <- rbind(cbind(prior$Pi$constant$V_i, matrix(0, n_Pi, n_B)), cbind(matrix(0, n_B, n_Pi), prior$B$constant[[2]]))
  if (tvp[1]) {
    Pi_Q_df_post <- t + data$priors$Pi$tvp[[1]]
    Pi_Q_V_prior <- data$priors$Pi$tvp[[2]]
  }
  if (tvp[2]) {
    B_Q_df_post <- t + data$priors$B$tvp[[1]]
    B_Q_V_prior <- data$priors$B$tvp[[2]]
  }
  if (tvp[1] & tvp[2]) {
    Pi_B_Q_df_post <- rbind(Pi_Q_df_post, B_Q_df_post)
    Pi_B_Q_V_prior <- rbind(cbind(Pi_Q_V_prior, matrix(0, n_Pi, n_B)), cbind(matrix(0, n_B, n_Pi), B_Q_V_prior))
  }
  if (shrinkage){
    shrinkage_type <- prior$Shrinkage$type
    if (shrinkage_type == "BVS"){
      Pi_B_restricted_variables <- prior$Shrinkage$spec$B[[1]] + n_Pi
      Pi_B_lpr_include_prior <- matrix(c(rep(1, n_Pi), prior$Shrinkage$spec$B[[2]]))
      Pi_B_lpr_exclude_prior <- matrix(c(rep(1, n_Pi), prior$Shrinkage$spec$B[[3]]))
    }
  }
  
  # Covariance matrix
  ## Constant covariance
  if (structural) {
    A0_mu_prior <- prior$A0$constant[[1]]
    A0_V_i_prior <- prior$A0$constant[[2]]
    #if (tvp[3]) {
    #  A0_Q_df_post <- t + prior$A0$tvp[[1]]
    #  A0_Q_V_prior <- prior$A0$tvp[[2]]
    #}
    #if (shrinkage){
    #  if (shrinkage_type == "BVS"){
    #    A0_lpr_include_prior <- prior$Shrinkage$spec$A0[[2]]
    #    A0_lpr_exclude_prior <- prior$Shrinkage$spec$A0[[3]]
    #  }
    # }
    
    if (sv) {
      stoch_vol <- prior$Sigma$sv
    }
    
    Sigma_df_post <- t + prior$Omega$constant[[1]]
    Sigma_V_prior <- prior$Omega$constant[[2]] 
    #}
  } else {
    Omega_df_post <- t + prior$Omega$constant[[1]]
    Omega_V_prior <- prior$Omega$constant[[2]]
  }
  
  #### Initialise ####
  pos_Pi <- 1:n_Pi
  pos_B <- n_Pi + 1:n_B
  ect_x <- rbind(ect, x)
  
  if (tvp[1] | tvp[2]) {
    y_0 <- y * NA
    y_dot <- y * NA
    Pi_B <- matrix(0, n_Pi_B, t + 1)
    Pi_B_filter <- Pi_B
  }
  if (tvp[1]) {
    Pi_Q <- diag(1, n_Pi)
    Pi_res <- cbind(matrix(diag(Pi_Q_V_prior)), Pi_Q_df_post - t)
    diag(Pi_Q) <- 1 / apply(Pi_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    Pi_T <- diag(.999, n_Pi)
  }
  if (tvp[2]) {
    B_Q <- diag(1, n_B)
    B_res <- cbind(matrix(diag(B_Q_V_prior)), B_Q_df_post - t)
    diag(B_Q) <- 1 / apply(B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    B_T <- diag(1, n_B)
  }
  if (tvp[1] & tvp[2]) {
    Pi_B_Q <- rbind(cbind(Pi_Q, matrix(0, n_Pi, n_B)), cbind(matrix(0, n_B, n_Pi), B_Q))
    Pi_B_T <- rbind(cbind(Pi_T, matrix(0, n_Pi, n_B)), cbind(matrix(0, n_B, n_Pi), B_T))
  }
  if (tvp[1] | tvp[2] | tvp_covar | shrinkage) {
    Z_Pi <- kronecker(t(ect), diag(1, n))
    Z_B <- kronecker(t(x), diag(1, n))
    Z_Pi_B <- cbind(Z_Pi, Z_B)
    if (shrinkage) {
      Z_Pi_B_raw <- Z_Pi_B
    }
  }
  if (shrinkage) {
    Pi_B_Gamma <- diag(1, n_Pi_B)
  }
  ## Variance block
  if (structural) {
    Sigma <- diag(apply(y, 1, stats::var), n)
    Sigma_i <- solve(Sigma)
    if (sv) {
     Sigma_temp <- Sigma
     Sigma_i_temp <- Sigma_i
     Sigma <- matrix(0, n * t, n)
     Sigma_i <- Sigma
     for (i in 1:t) {
       Sigma[(i - 1) * n + 1:n,] <- Sigma_temp
       Sigma_i[(i - 1) * n + 1:n,] <- Sigma_i_temp
     }
     rm(list = c("Sigma_temp", "Sigma_i_temp"))
    }
    
    #if (tvp[3]) {
    #   A0 <- t(matrix(diag(1, n), n, n * (t + 1)))
    #   A0_filter <- matrix(0, n_A0, t + 1)
    #   A0_Q <- diag(.0001, n_A0)
    # } else {
    A0 <- diag(1, n)
    #}
    #   if (shrinkage) {
    #     if (structural_est) {
    #       Gamma_A0 <- diag(1, n_A0)
    #     }
    #   }
    #}
    Omega <- Sigma
    Omega_i <- Sigma_i
  } else {
    Omega <- diag(apply(y, 1, stats::var), n)
    Omega_i <- solve(Omega)
  }
  
  #### Data containers ####
  # Cointegration parameters
  if (tvp[1]) {
    draws_Pi <- array(NA, c(n * n_ect, t, store))
  } else {
    draws_Pi <- array(NA, c(n * n_ect, 1, store))
  }
  # Non-cointegration parameters
  if (tvp[2]) {
    draws_B <- array(NA, c(n_B, t, store))
  } else {
    draws_B <- array(NA, c(n_B, 1, store))
  }
  if (shrinkage) {
    draws_B_shrinkage <- matrix(NA, n_B, iterations - burnin)
    #if (structural_est) {
    #  draws_A0_shrinkage <- matrix(NA, n_A0, store)
    #}
  }
  # Covariance matrix
  if (structural) {
    if (structural_est) {
      if (tvp[3]) {
        draws_A0 <- array(NA, c(n * n, t, store))
      } else {
        draws_A0 <- array(NA, c(n * n, 1, store))
      }
    }
    if (sv) {
      draws_Sigma <- array(NA, c(n * n, t, store))
    } else {
      draws_Sigma <- array(NA, c(n * n, 1, store))
    }
  }
  if (tvp_covar) {
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
    # Draw parameters
    if (shrinkage) {
      Z_Pi_B <- Z_Pi_B_raw %*% Pi_B_Gamma
    }
    if (tvp[1] | tvp[2]) {
      for (i in 1:t) {
        y_0[, i] <- y[, i] - Z_Pi_B[(i - 1) * n + 1:n,] %*% Pi_B_filter[, i]
      }
      Pi_B_0 <- posterior_normal_sur(y_0, Z_Pi_B, Omega_i, Pi_B_mu_prior, Pi_B_V_i_prior)
      y_filter <- y - matrix(Z_Pi_B %*% Pi_B_0, n)
      
      if (tvp[1] & tvp[2]) {
        stop("No!")
        Pi_B_filter <- dk(y_filter, Z_Pi_B, Omega, Pi_B_Q, Pi_B_T, matrix(0, n_Pi_B), Pi_B_Q)
        Pi_B <- Pi_B_filter + matrix(Pi_B_0, n_Pi_B, t + 1)
        Pi <- Pi_B[pos_Pi,]
        B <- Pi_B[pos_B,]
        
        Pi_B_res <- cbind(matrix(diag(Pi_B_Q_V_prior) + rowSums((Pi_B[,-1] - Pi_B[, -(t + 1)])^2)), Pi_B_Q_df_post)
        diag(Pi_B_Q) <- 1 / apply(Pi_B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
      } else {
        if (tvp[2]) {
          B_filter <- dk(y_filter, Z_Pi_B[, pos_B], Omega, B_Q, B_T, matrix(0, n_B), B_Q)
          Pi_B_filter[pos_B,] <- B_filter
          
          Pi_B <- Pi_B_filter + matrix(Pi_B_0, n_Pi_B, t + 1)
          
          if (shrinkage) {
            Pi_B_Gamma <- bvs(y, Z_Pi_B_raw, Pi_B[,-(t + 1)], Pi_B_Gamma, Omega_i, Pi_B_restricted_variables,
                              Pi_B_lpr_include_prior, Pi_B_lpr_exclude_prior)
            Pi_B <- Pi_B_Gamma %*% Pi_B 
          }
          
          Pi <- Pi_B_0[pos_Pi, 1]
          B <- Pi_B[pos_B,]
          
          B_res <- cbind(matrix(diag(B_Q_V_prior) + rowSums((B[,-1] - B[, -(t + 1)])^2)), B_Q_df_post)
          diag(B_Q) <- 1 / apply(B_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
        }
      }
      for (i in 1:t) {
        y_dot[, i] <- y[, i] - matrix(Pi_B[, i], n) %*% ect_x[,i]
      }
    } else {
      if (tvp_covar | shrinkage) {
        Pi_B <- posterior_normal_sur(y, Z_Pi_B, Omega_i, Pi_B_mu_prior, Pi_B_V_i_prior)
        
        if (shrinkage) {
          Pi_B_Gamma <- bvs(y, Z_Pi_B_raw, Pi_B, Pi_B_Gamma, Omega_i, Pi_B_restricted_variables,
                            Pi_B_lpr_include_prior, Pi_B_lpr_exclude_prior)
          Pi_B <- Pi_B_Gamma %*% Pi_B 
        }
      } else {
        Pi_B <- posterior_normal(y, ect_x, Omega_i, Pi_B_mu_prior, Pi_B_V_i_prior)
      }
      y_dot <- y - matrix(Pi_B, n) %*% ect_x
      Pi <- Pi_B[pos_Pi,]
      B <- Pi_B[pos_B,]
    }

    # Structural parameters
    if (structural) {
      if (structural_est) {
        A0_temp <- draw_structural_parameters(y_dot, Sigma_i, A0_mu_prior, A0_V_i_prior, tvp[3])
        A0 <- A0_temp[[1]]
        y_dot_dot <- A0_temp[[2]]
      }
      
      # Volatility states
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
          Sigma[i + 0:(t - 1) * n, i] <- h
          Sigma_i[i + 0:(t - 1) * n, i] <- 1 / h
        }
      } else {
        Sigma_res <- cbind(diag(Sigma_V_prior) + rowSums(y_dot_dot^2), Sigma_df_post)
        diag(Sigma_i) <- apply(Sigma_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
        diag(Sigma) <- 1 / diag(Sigma_i)
      }
      
      # Combine structural parameters with volatilities
      if (tvp[3]) {
        for (i in 1:t) {
         A0_i_temp <- solve(A0[(i - 1) * n + 1:n, ])
         if (sv) {
           Omega[(i - 1) * n + 1:n,] <- A0_i_temp %*% tcrossprod(Sigma[(i - 1) * n + 1:n,], A0_i_temp)
         } else {
           Omega[(i - 1) * n + 1:n,] <- A0_i %*% tcrossprod(Sigma, A0_i)
         }
         Omega_i[(i - 1) * n + 1:n,] <- solve(Omega[(i - 1) * n + 1:n,])
        }
      } else {
        A0_i <- solve(A0)
        if (sv) {
          for (i in 1:t) {
            Omega[(i - 1) * n + 1:n, ] <- A0_i %*% tcrossprod(Sigma[(i - 1) * n + 1:n, ], A0_i)
            Omega_i[(i - 1) * n + 1:n, ] <- solve(Omega[(i - 1) * n + 1:n, ]) 
          }
        } else {
          Omega <- A0_i %*% tcrossprod(Sigma, A0_i)
          Omega_i <- solve(Omega) 
        }
      }
    } else {
      Omega_i <- rWishart(1, Omega_df_post, solve(tcrossprod(y_dot)))[,, 1] 
      Omega <- solve(Omega_i)
    }
    
    #### Chunks ####

    # Korobilis (2013)
    #   if (shrinkage) {
    #     B_Gamma <- bvs(y_tilde, Z_x, B[, -(t + 1)], B_Gamma, Omega_i, B_restricted_variables, B_lpr_include_prior, B_lpr_exclude_prior)
    #     B <- B_Gamma %*% B
    #   }
    #       # Initialise
    #       for (i in 1:t) {
    #         y_hat_0[, i] <- y_hat[, i] - Z_Pi[(i - 1) * n + 1:n,] %*% Pi_filter[, i]
    #       }
    #       Pi_0 <- posterior_normal_sur(y_hat_0, Z_Pi, Omega_i, Pi_mu_prior, Pi_V_i_prior)
    #       y_hat_filter <- y_hat - matrix(Z_Pi %*% Pi_0, n)
    #       
    #       # Durbin & Koopman (2002)
    #       Pi_filter <- dk(y_hat_filter, Z_Pi, Omega, Pi_Q, Pi_T, Pi_a0, Pi_Q)
    #       Pi_temp <- Pi_filter + matrix(Pi_0, n_alpha, t + 1)
    #       for (i in 1:t) {
    #         Pi[(i - 1) * n + 1:n,] <- matrix(Pi_temp[,i], n)
    #       }
    #       
    #       # Covariance matrix of the coefficients
    #       Pi_res <- cbind(matrix(diag(Pi_Q_V_prior) + rowSums((Pi_temp[,-1] - Pi_temp[, -(t + 1)])^2)), Pi_Q_df_post)
    #       diag(Pi_Q) <- 1 / apply(Pi_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    #     } 
    #     for (i in 1:t) {
    #       y_dot[, i] <- y_hat[, i] - matrix(Pi[, i], n) %*% ect[, i]
    #     }
    # 
    # #### Error covariance matrix ####
    # if (structural) {
    #   #### Structural ####

    #       if (tvp[3]) {
    #         y_dot_temp_0 <- y_dot_temp * 0
    #         A0_filter_temp <- matrix(A0_filter[p1,], i - 1)
    #         if (sv) {
    #           Sigma2_temp <- matrix(Sigma2[(0:(t - 1) * n) + i,i])
    #         } else {
    #           Sigma2_temp <- matrix(diag(Sigma2)[i]) 
    #         }
    #         
    #         # Initialise
    #         for (j in 1:t) {
    #           y_dot_temp_0[, j] <- y_dot_temp[, j] - x_dot_temp[, j] %*% A0_filter_temp[, j]
    #         }
    #         A0_0 <- matrix(posterior_normal_sur(y_dot_temp_0, t(x_dot_temp), Sigma2_i_temp, A0_mu_prior_temp, A0_V_i_prior_temp), 1)
    #         y_dot_temp_filter <- y_dot_temp - A0_0 %*% x_dot_temp
    #         
    #         # Durbin & Koopman (2002)
    #         A0_filter_temp <- dk(y_dot_temp_filter, t(x_dot_temp), Sigma2_temp, matrix(A0_Q[p1, p1], i - 1), diag(1, i - 1),
    #                              matrix(0, i - 1), matrix(A0_Q[p1, p1], i - 1))
    #         A0_temp <- A0_filter_temp + matrix(A0_0, i - 1, t + 1)
    #         
    #         # Korobilis (2013)
    #         if (shrinkage) {
    #           Gamma_A0[p1, p1] <- bvs(y_dot_temp, t(x_dot_temp), matrix(A0_temp[, -(t + 1)], i - 1), matrix(Gamma_A0[p1, p1], i - 1),
    #                                   Sigma2_i_temp, matrix(1:(i - 1)), A0_lpr_include_prior, A0_lpr_exclude_prior)
    #           A0_temp <- Gamma_A0[p1, p1] %*% A0_temp
    #         }
    #         
    #         # Covariance matrix of the coefficients
    #         A0_res <- cbind(matrix(diag(A0_Q_V_prior)[p1] + rowSums(matrix(A0_temp[,-1] - A0_temp[, -(t + 1)], i - 1)^2)), A0_Q_df_post[p1])
    #         diag(A0_Q)[p1] <- 1 / apply(A0_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
    #         
    #         A0[i + (0:t * n), 1:(i - 1)] <- t(A0_temp)
    #       } else {

    #     
    #     #### Volatility state ####
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
    #   }
    # 
    #### Post-burnin procedures ####
    if (draw > burnin){
      pos_draw <- draw - burnin
      
      if (tvp[1]) {
        draws_Pi[,, pos_draw] <- Pi[, -(t + 1)]
      } else {
        draws_Pi[, 1, pos_draw] <- Pi
      }
      
      if (tvp[2]) {
        draws_B[, , pos_draw] <- B[, -(t + 1)]
      } else {
        draws_B[, 1, pos_draw] <- B
      }
      
      if (shrinkage) {
        draws_B_shrinkage[, pos_draw] <- diag(Pi_B_Gamma)[pos_B]
        #if (structural_est) {
        #  draws_A0_shrinkage[, pos_draw] <- diag(Gamma_A0)
        #}
      }
      
      if (structural) {
        if (structural_est) {
          #if (tvp[3]) {
          #  for (i in 1:t) {
          #    draws_A0[, i, pos_draw] <- A0[, -(t + 1)]
          #  }
          #} else {
            draws_A0[, 1, pos_draw] <- A0
          #}
        }
        if (sv) {
          for (i in 1:t) {
            draws_Sigma[, i, pos_draw] <- Sigma[(i - 1) * n + 1:n,]
          }
        } else {
          draws_Sigma[, 1, pos_draw] <- Sigma
        }
      }
      
      if (tvp_covar) {
        for (i in 1:t) {
          draws_Omega[, i, pos_draw] <- Omega[(i - 1) * n + 1:n,] 
        }
      } else {
        draws_Omega[, 1, pos_draw] <- Omega
      }
      
      draws_LL[, pos_draw] <- getLL(y_dot, Omega, Omega_i)
      ###########################################################################
      if (FALSE) {
        #### Print ####
        p <- estimation_data$domestic["lag"]
        graphics::par(mfcol = c(2, 2))
        # 
        stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter y")
        
        stats::plot.ts(cbind(0, draws_B[n * n * p + 1, 1, ]), plot.type = "single", col = "black", ylab = "Contemp y")
        # 
        stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                       plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic y")
        stats::plot.ts(cbind(0, draws_B[n_B - 2, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic y")
        # 
        #stats::plot.ts(cbind(0, draws_Pi[n_ect - 2, 1, ]), plot.type = "single", col = "black", ylab = "Pi y")
        # 
        #stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + n + 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter Dp")
        #stats::plot.ts(cbind(0, draws_B[n * n * p + n + 2, 1, ]), plot.type = "single", col = "black", ylab = "Contemp Dp")
        # 
        #stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 1, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic Dp")
        #stats::plot.ts(cbind(0, draws_B[n_B - 1, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic Dp")
        # 
        #stats::plot.ts(cbind(0, draws_Pi[n_ect - 1, 1, ]), plot.type = "single", col = "black", ylab = "Pi Dp")
        # 
        #stats::plot.ts(cbind(0, t(apply(draws_B[n * n * p + 2 * n + 3, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter r")
        #stats::plot.ts(cbind(0, draws_B[n * n * p + 2 * n + 3, 1, ]), plot.type = "single", col = "black", ylab = "Contemp r")
        # 
        #stats::plot.ts(cbind(0, t(apply(draws_B[n_B, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Deterministic r")
        #stats::plot.ts(cbind(0, draws_B[n_B, 1, ]), plot.type = "single", col = "black", ylab = "Deterministic r")
        # 
        #stats::plot.ts(cbind(0, draws_Pi[n_ect, 1, ]), plot.type = "single", col = "black", ylab = "Pi r")
        #stats::plot.ts(cbind(0, t(apply(draws_B[n_B - 2, , (ceiling(pos_draw) * .5):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #               plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter")
        # 
        # 
        # stats::plot.ts(log(t(apply(draws_Sigma[1, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #         plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Log-Volatility y")
        # 
        # stats::plot.ts(log(t(apply(draws_Sigma[5, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #                plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Log-Volatility Dp")
        # 
        # stats::plot.ts(log(t(apply(draws_Sigma[9, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #                plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Log-Volatility r")
        # 
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
  
  tw <- sample(x = 1:store, size = store / thin, replace = FALSE)
  tw <- tw[order(tw)]
  
  draws <- length(tw)
  A0_t <- ifelse(tvp[3], t, 1)
  if (structural_est) {
    if (tvp[3]) {
      A0 <- draws_A0[,, tw]
    } else {
      A0 <- array(draws_A0[, 1, tw], dim = c(n * n, 1, draws))
    }
  } else {
    A0 <- array(diag(1, n), dim = c(n * n, A0_t, draws))
  }
  dimnames(A0) <- list(paste(names_y, rep(names_y, each = n), sep = "_"), NULL, NULL)
  
  B_t <- ifelse(tvp[2], t, 1)
  if (tvp[2]) {
    B <- draws_B[,, tw]
  } else {
    B <- array(draws_B[, 1, tw], dim = c(n_B, 1, draws))
  }
  
  Pi_t <- ifelse(tvp[1], t, 1)
  if (tvp[1]) {
    Pi <- draws_Pi[,, tw]
  } else {
    Pi <- array(draws_Pi[,, tw], dim = c(n * n_ect, 1, draws)) 
  }
  
  if (structural) {
    if (sv) {
      Sigma <- draws_Sigma[,, tw]   
    } else {
      Sigma <- array(draws_Sigma[, 1, tw], dim = c(n * n, 1, draws))
    }
    dimnames(Sigma) <- list(paste(names_y, rep(names_y, each = n), sep = "_"), NULL, NULL)
  }
  
  if (tvp_covar) {
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
  test_stat <- c("AIC" = 2 * sum(log(rowMeans(draws_LL))) - (n_Pi_B) * 2,
                 "BIC" = 2 * sum(log(rowMeans(draws_LL))) - (n_Pi_B) * log(t),
                 "HQ" = 2 * sum(log(rowMeans(draws_LL))) - (n_Pi_B) * 2 * log(log(t)))
  
  ##### Prepare output object ####
  data$specs$tvp <- c((tvp[1] | tvp[2]), tvp[3])
  result <- c(data, list(coefs = list(A0 = A0, A_d = A_d, A_s_0 = A_s_0, A_s = A_s, A_g = A_g, A_det = det, Omega = Omega),
                         criteria = list("Criteria" = test_stat)))
  
  if (structural) {
    if (sv) {
      result$coefs$Sigma <- draws_Sigma[,,tw] 
    } else {
      result$coefs$Sigma <- array(draws_Sigma[,,tw], c(n * n, 1, draws))
    }
  }
  
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
