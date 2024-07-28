
# Algorithm for Bayesian inference of VEC models with time varying parameters

.bvectvpalg <- function(object) {
  
  # Data ----
  y <- t(object[["data"]][["y"]])
  yvec <- matrix(y)
  w <- object[["data"]][["w"]]
  z <- object[["data"]][["z"]]
  
  k <- nrow(y)
  r <- object[["model"]][["rank"]]
  tt <- ncol(y)
  if (is.null(z)) {
    n_a <- 0
  } else {
    n_a <- ncol(z)
  }
  use_a <- n_a > 0
  
  vec_tt <- matrix(1, tt)
  diag_tt <- Matrix::Diagonal(tt, x = 1)
  
  n_alpha <- r * k
  use_nonalpha <- n_alpha < n_a
  k_ect <- ncol(w)
  n_beta <- r * k_ect
  use_rr <- r > 0
  
  # Priors and initial values ----
  
  ## Non-cointegration ----
  use_a_bvs <- FALSE
  if (use_a) {
    zz_a <- sur_const_to_tvp(z, k, tt)
    rho <- object[["priors"]][["cointegration"]][["rho"]]
    a_v_i <- object[["initial"]][["a_v_i"]]
    a_v_prior_shape <- object[["priors"]][["a"]][["shape"]]
    a_v_prior_rate <- object[["priors"]][["a"]][["rate"]]
    a_init <- object[["initial"]][["a_init"]]
    a_init_prior_mu <- object[["priors"]][["a"]][["mu"]]
    a_init_prior_vi <- object[["priors"]][["a"]][["v_i"]]
    hh_a <- create_first_difference_matrix(n_a, tt)
    if (use_nonalpha) {
      pos_non_alpha <- (n_alpha + 1):n_a
      pos_non_alpha_large <- rep(pos_non_alpha, tt) + rep(0:(tt - 1), each = length(pos_non_alpha)) * n_a 
    }
    
    if (!is.null(object[["model"]][["varselect"]])) {
      if (object[["model"]][["varselect"]] == "BVS" & !is.null(object[["priors"]][["a"]][["bvs"]])) {
        use_a_bvs <- TRUE
        a_bvs_prior <- object[["priors"]][["a"]][["bvs"]][["inprior"]]
        a_bvs_include <- object[["priors"]][["a"]][["bvs"]][["include"]]
        a_bvs_lambda <- Matrix::Diagonal(n_a, 1)
      }
    }
    
    ## Cointegration ----
    if (use_rr) {
      beta <- object[["initial"]][["beta"]]
      beta_init <- object[["initial"]][["beta_init"]]
      z_b <- matrix(NA, k * tt, n_beta)
      hh_b <- create_first_difference_matrix(n_beta, tt)
      beta_v_i <- Matrix::Diagonal(n = n_beta, x = 1)
      hb_vi_hb <- crossprod(hh_b, kronecker(diag_tt, beta_v_i)) %*% hh_b
      beta_init_prior_mu <- object[["priors"]][["cointegration"]][["mu"]]
      beta_init_prior_vi <- object[["priors"]][["cointegration"]][["v_i"]]
      # Since those matrices are fixed they can be defined here
      beta_init_post_v <-  beta_init_prior_vi + beta_v_i
    }
  }
  
  ## Covariances ----
  use_psi <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  use_psi_bvs <- FALSE
  if (use_psi) {
    n_psi <- k * (k - 1) / 2
    psi_v_i <- object[["initial"]][["psi_v_i"]]
    psi_v_prior_shape <- object[["priors"]][["psi"]][["shape"]]
    psi_v_prior_rate <- object[["priors"]][["psi"]][["rate"]]
    psi_init <- object[["initial"]][["psi_init"]]
    psi_init_prior_mu <- object[["priors"]][["psi"]][["mu"]]
    psi_init_prior_vi <- object[["priors"]][["psi"]][["v_i"]]
    hh_psi <- create_first_difference_matrix(n_psi, tt)
    
    psi_pos_x <- 1:(k * tt)
    psi_pos_y <- 1:(k * tt)
    for (j in 1:tt) {
      for (i in 1:(k - 1)) {
        psi_pos_x <- append(psi_pos_x, 1:i + (j - 1) * k)
      } 
      for (i in 2:k) {
        psi_pos_y <- append(psi_pos_y, rep(i, i - 1) + (j - 1) * k)
      } 
    }
    Psi <- Matrix::spMatrix(nrow = k * tt,
                            ncol = k * tt,
                            i = psi_pos_y,
                            j = psi_pos_x,
                            x = c(rep(1, k * tt), rep(2, n_psi * tt)))
    psi_pos <- Matrix::which(Psi == 2)
    Psi[psi_pos] <- 0
    
    if (!is.null(object[["model"]][["varselect"]])) {
      if (object[["model"]][["varselect"]] == "BVS" & !is.null(object[["priors"]][["psi"]][["bvs"]])) {
        use_psi_bvs <- TRUE
        psi_bvs_prior <- object[["priors"]][["psi"]][["bvs"]][["inprior"]]
        psi_bvs_include <- object[["priors"]][["psi"]][["bvs"]][["include"]]
        psi_bvs_lambda <- Matrix::Diagonal(n_psi, 1)
      }
    }
  }
  
  ## Error term ----
  u <- yvec
  use_sv <- object[["model"]][["error"]] %in% c("sv", "sv+covar")
  use_gamma <- object[["model"]][["error"]] %in% c("gamma", "gamma+covar")
  if (use_sv) {
    h_init_prior_mu <- object[["priors"]][["sigma"]][["mu"]]
    h_init_prior_vi <- object[["priors"]][["sigma"]][["v_i"]]
    h_v_prior_shape <- object[["priors"]][["sigma"]][["shape"]]
    h_v_prior_rate <- object[["priors"]][["sigma"]][["rate"]]
    h <- object[["initial"]][["h"]]
    h_init <- object[["initial"]][["h_init"]]
    h_v <- object[["initial"]][["h_state_variance"]]
    h_offset <- object[["initial"]][["h_offset"]]
    omega_i_diag <- Matrix::Diagonal(k * tt, 1)
    Matrix::diag(omega_i_diag) <- 1 / exp(matrix(t(h)))
    sigma_i_diag <- omega_i_diag
  } else {
    if (use_gamma) {
      omega_i <- object[["initial"]][["sigma_i"]]
      omega_i_diag <- kronecker(diag_tt, omega_i)
      omega_i_prior_shape <- object[["priors"]][["sigma"]][["shape"]]
      omega_i_prior_rate <- object[["priors"]][["sigma"]][["rate"]]
    } else {
      sigma_prior_scale <- object[["priors"]][["sigma"]][["scale"]]
      sigma_post_df <- object[["priors"]][["sigma"]][["df"]] + tt 
    }
  }
  if (use_sv | use_psi) {
    # Get positions of the elements of a block diagonal matrix
    temp <- list()
    for (i in 1:tt) {
      temp[[i]] <- matrix(1, k, k)
    }
    temp <- Matrix::bdiag(temp)
    pos_sigma <- Matrix::which(temp == 1)
    rm(temp)
  }
  if (!use_sv) {
    sigma_i <- object[["initial"]][["sigma_i"]]
    sigma_i_diag <- kronecker(diag_tt, sigma_i) 
  }
  
  
  # Storage objects ----
  iterations <- object[["model"]][["iterations"]]
  burnin <- object[["model"]][["burnin"]]
  draws <- iterations + burnin
  
  draws_a <- NULL
  draws_a_sigma <- NULL
  draws_beta <- NULL
  if (use_a) {
    draws_a <- matrix(NA, iterations, n_a * tt)
    draws_a_sigma <- matrix(NA, iterations, n_a)
    if (use_a_bvs) {
      draws_a_lambda <- matrix(NA, iterations, n_a)
    }
    draws_beta <- matrix(NA, iterations, n_beta * tt)
  }
  
  if (use_psi | use_sv) {
    draws_sigma <- matrix(NA, iterations, k * k * tt)
  } else {
    draws_sigma <- matrix(NA, iterations, k * k) 
  }
  
  if (use_psi_bvs) {
    draws_psi_lambda <- matrix(NA, iterations, n_psi)
  }
  
  pb <- utils::txtProgressBar(style = 3)
  
  # Gibbs sampler ----
  for (draw in 1:draws) {
    
    if (use_a) {
      
      ## Draw non-cointegration parameters ----
      
      ### Draw states ----
      h_vi_h <- t(hh_a) %*% kronecker(diag_tt, a_v_i) %*% hh_a
      a_v_post <- h_vi_h + t(zz_a) %*% sigma_i_diag %*% zz_a
      a_mu_post <- matrix(solve(a_v_post, h_vi_h %*% kronecker(vec_tt, a_init) + t(zz_a) %*% sigma_i_diag %*% yvec))
      a <- matrix(a_mu_post + solve(chol(a_v_post), rnorm(n_a * tt)))
      
      ### Variable selection ----
      if (use_a_bvs) {
        a_bvs_lambda <- post_bvs(yvec, z, a, k, n_a, a_bvs_lambda, sigma_i_diag, a_bvs_prior, a_bvs_include) 
        a <- matrix(kronecker(diag_tt, a_bvs_lambda) %*% a)
        a_init <- matrix(a_bvs_lambda %*% a_init)
      }
      
      ### Draw state variances ----
      a_v_i <- post_gamma_state_variance(a = a, a_init = a_init,
                                         shape_prior = a_v_prior_shape,
                                         rate_prior = a_v_prior_rate,
                                         inverse = TRUE) 
      
      ### Draw initial values ----
      a_init_post_v <- a_init_prior_vi + a_v_i
      a_init_post_mu <- solve(a_init_post_v, a_init_prior_vi %*% a_init_prior_mu + a_v_i %*% a[1:n_a])
      a_init <- matrix(a_init_post_mu + solve(chol(a_init_post_v), rnorm(n_a)))
      
      ## Draw cointegration parameters ----
      if (use_rr) {
        
        ### Prepare data ----
        if (use_nonalpha) {
          ystar <- yvec - zz_a[, pos_non_alpha_large] %*% a[pos_non_alpha_large] 
        } else {
          ystar <- yvec
        }
        
        for (i in 1:tt) {
          z_b[(i - 1) * k + 1:k,] <- kronecker(matrix(a[(i - 1) * n_a + 1:n_alpha], k), matrix(w[i, ], 1))
        }
        zz_b <- sur_const_to_tvp(z_b, k, tt)
        
        ### Draw states ----
        beta_v_post <- hb_vi_hb + t(zz_b) %*% sigma_i_diag %*% zz_b
        beta_mu_post <- solve(beta_v_post, hb_vi_hb %*% kronecker(vec_tt, beta_init) + t(zz_b) %*% sigma_i_diag %*% ystar)
        beta <- matrix(beta_mu_post + solve(chol(beta_v_post), rnorm(n_beta * tt)))
        
        ### Draw initial state ----
        beta_init_post_mu <- matrix(solve(beta_init_post_v, beta_init_prior_vi %*% beta_init_prior_mu + beta_v_i %*% beta[1:n_beta]))
        beta_init <- beta_init_post_mu + solve(chol(beta_init_post_v), rnorm(n_beta))
        
        ### Update z data table
        for (i in 1:tt) {
          z[(i - 1) * k + 1:k, 1:n_alpha] <- kronecker(t(crossprod(matrix(beta[(i - 1) * n_beta + 1:n_beta], k_ect), w[i, ])), diag(1, k))
        }
        zz_a <- sur_const_to_tvp(z, k, tt)
      }
      
      u <- matrix(yvec - zz_a %*% a)
    }
    
    
    ## Draw covariance parameters ----
    if (use_psi) {
      
      psi <- post_normal_covar_tvp(u, omega_i_diag, k, psi_v_i, psi_init)
      
      ### Variable selection ----
      if (use_psi_bvs) {
        
        psi_temp <- covar_prepare_data(y = u, omega_i = omega_i_diag, k = k, tt = tt, tvp = FALSE)
        
        psi_bvs_lambda <- post_bvs(y = psi_temp$y,
                                   z = as.matrix(psi_temp$z),
                                   a = psi,
                                   k = k - 1,
                                   m = n_psi,
                                   lambda = psi_bvs_lambda,
                                   sigma_i = psi_temp$omega_i,
                                   prob_prior = psi_bvs_prior,
                                   include = psi_bvs_include) 
        
        psi <- matrix(kronecker(diag_tt, psi_bvs_lambda) %*% psi)
        psi_init <- matrix(psi_bvs_lambda %*% psi_init)
      }
      
      ### Draw state variances ----
      psi_v_i <- post_gamma_state_variance(a = psi, a_init = psi_init,
                                           shape_prior = psi_v_prior_shape,
                                           rate_prior = psi_v_prior_rate,
                                           inverse = TRUE) 
      
      ### Draw initial values ----
      psi_init_post_v <- psi_init_prior_vi + psi_v_i
      psi_init_post_mu <- solve(psi_init_post_v, psi_init_prior_vi %*% psi_init_prior_mu + psi_v_i %*% psi[1:n_psi])
      psi_init <- matrix(psi_init_post_mu + solve(chol(psi_init_post_v), rnorm(n_psi)))
      
      Psi[psi_pos] <- psi # Turn psi vector into block diagonal matrix
      u <- matrix(Psi %*% u) # Calculate errors
    }
    
    ## Draw error variances ----
    if (use_sv) {
      
      ### Draw states ----
      h <- stochvol_ksc1998(y = t(matrix(u, k, tt)),
                            h = h,
                            sigma = h_v,
                            h_init = h_init,
                            constant = h_offset)
      
      ### Draw state variances ----
      h_v_i <- post_gamma_state_variance(a = matrix(t(h)),
                                         a_init = h_init,
                                         shape_prior = h_v_prior_shape,
                                         rate_prior = h_v_prior_rate,
                                         inverse = TRUE)
      
      h_v <- 1 / matrix(Matrix::diag(h_v_i))
      
      ### Draw initial values ----
      h_init_post_v <- h_init_prior_vi + h_v_i
      h_init_post_mu <- solve(h_init_post_v, h_init_prior_vi %*% h_init_prior_mu + h_v_i %*% h[1,])
      h_init <- matrix(h_init_post_mu + solve(chol(h_init_post_v), rnorm(k)))
      
      
      Matrix::diag(omega_i_diag) <- 1 / exp(matrix(t(h)))
      if (use_psi) {
        sigma_i_diag <- t(Psi) %*% omega_i_diag %*% Psi
      } else {
        sigma_i_diag <- omega_i_diag;
      }
      
    } else {
      if (use_gamma) {
        
        omega_i <- post_gamma_measurement_variance(u = u,
                                                   shape_prior = omega_i_prior_shape,
                                                   rate_prior = omega_i_prior_rate,
                                                   inverse = TRUE)
        
        if (use_psi) {
          omega_i_diag <- kronecker(diag_tt, omega_i)
          sigma_i_diag <- t(Psi) %*% omega_i_diag %*% Psi
        } else {
          sigma_i <- omega_i
          sigma_i_diag <- kronecker(diag_tt, sigma_i)
        }
      } else {
        sigma_i <- stats::rWishart(1, sigma_post_df, solve(sigma_prior_scale + tcrossprod(matrix(u, k))))[,,1]
        sigma_i_diag <- kronecker(diag_tt, sigma_i) 
      }
    }
    
    ## Save draws
    if (draw > burnin) {
      
      pos_draw <- draw - burnin
      
      if (use_a) {
        draws_a[pos_draw, ] <- a 
        draws_a_sigma[pos_draw, ] <- 1 / Matrix::diag(a_v_i)
        if (use_a_bvs) {
          draws_a_lambda[pos_draw, ] <- Matrix::diag(a_bvs_lambda)
        }
        
        if (use_rr) {
          draws_beta[pos_draw,] <- beta
        } 
      }
      
      if (use_psi | use_sv) {
        draws_sigma[pos_draw, ] <- solve(sigma_i_diag)[pos_sigma]
        if (use_psi_bvs) {
          draws_psi_lambda[pos_draw, ] <- Matrix::diag(psi_bvs_lambda)
        }
      } else {
        draws_sigma[pos_draw, ] <- matrix(solve(sigma_i)) 
      }
    }
    
    if (draw %% 100) {
      utils::setTxtProgressBar(pb, draw / draws) 
    }
    
  }
  
  object <- c(object, list(posteriors = list(a = list(coeffs = draws_a,
                                                      sigma = draws_a_sigma),
                                             beta = list(coeffs = draws_beta),
                                             sigma = list(coeffs = draws_sigma))))
  
  if (use_a_bvs) {
    object[["posteriors"]][["a"]][["lambda"]] <- draws_a_lambda
  }
  if (use_psi_bvs) {
    object[["posteriors"]][["sigma"]][["lambda"]] <- draws_psi_lambda
  }
  
  return(object)
}