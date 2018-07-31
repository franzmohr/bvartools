#' Solve a GVAR Model
#' 
#' Combines the country model results to a global VAR model and solves it.
#' 
#' @param data a list containing the results of the country estimates.
#' @param ic a character specifying the information criterion used for model selection. Available options are "AIC", "BIC" (default) and "HQ".
# #' @param t an optional numeric specifying the period for which the GVAR model should be solved in case the models
#' were estimated with time varying parameters.
#' 
#' @details
#' If multiple country models were estimated with \code{\link{gvar_fit}}, the function will choose the model, which
#' maximises the specified information criterion specified by \code{ic}.
#' 
#' @return A list containing two components:
#' \item{Country.models}{a list of the results of all estimated country models as produced by \code{\link{estimate_country_model}}.}
#' \item{GVAR}{a list containing the data of the global model, the country-specific link matrices and a variable index as produced by the function \code{\link{country_models}}.}
#' 
#' @export
solve_gvar <- function(data, ic = "BIC", t = NULL){
  global.data <- data$global.data
  data <- data$country.data
  
  # Model selection
  if (!ic%in%c("AIC", "BIC", "HQ")) {
    stop("Information criterion not available. Please, select among 'AIC', 'BIC' or 'HQ'.")
  }
  W <- global.data$W
  index <- global.data$index
  
  r <- unlist(lapply(data, function(x){return(x$specs$rank$rank)}))
  Criteria <- t(matrix(unlist(lapply(data, function(x){
    if (is.null(x$coefs)) {g <- rep(NA, 3)} else {g <- as.numeric(x$criteria$Criteria)}
    return(g)
  })), 3))
  dimnames(Criteria)[[2]] <- c("AIC", "BIC", "HQ")
  criteria <- data.frame("Country" = names(r), "Rank" = r, Criteria, "Maximum" = FALSE)
  if (ic %in% c("AIC", "BIC", "HQ")){
    for (i in unique(names(data))) {
      c.pos <- which(criteria[, "Country"] == i)
      if (length(c.pos) > 1) {
        #pos.temp <- which(diff(criteria[c.pos, ic]) < 0)[1]
        pos.temp <- which.max(criteria[c.pos, ic])
        if (is.na(pos.temp)) {
          pos.temp <- which.max(criteria[c.pos, ic]) 
        }
        criteria[c.pos[pos.temp], "Maximum"] <- TRUE 
      } else {
        criteria[c.pos, "Maximum"] <- TRUE 
      }
    }
  } else {
    criteria[, "Maximum"] <- TRUE 
  }
  row.names(criteria) <- NULL
  pos.final <- which(criteria[, "Maximum"])
  c.final <- NULL
  n.final <- NULL
  for (i in pos.final){
    c.final <- c(c.final, list(data[[i]]))
    n.final <- c(n.final, names(data)[i])
  }
  names(c.final) <- n.final
  
  original.data <- data
  data <- c.final
  rm(c.final)
  
  # Solve GVAR model
  k <- dim(index)[1]
  n.c <- length(data)
  k.i <- unlist(lapply(data, function(x){length(x$specs$domestic.variables)}))
  k.s.i <- unlist(lapply(data, function(x){length(x$specs$foreign.variables)}))
  p.i <- unlist(lapply(data, function(x){x$specs$lags$lags$domestic}))
  p.star.i <- unlist(lapply(data, function(x){x$specs$lags$lags$foreign}))
  n.g <- unlist(lapply(data, function(x){g <- x$specs$global.variables; if (is.na(g)){return(0)} else {return(length(g))}}))
  p.global.i <- unlist(lapply(data, function(x){x$specs$lags$lags$global}))
  
  p <- max(p.i, p.star.i)
  if (all(is.na(p.global.i))){
    p.global <- NA
  } else {
    p.global <- max(stats::na.omit(p.global.i))
  }
  global <- !is.na(p.global)
  
  pos.i <- c(0, cumsum(k.i))
  
  draws <- dim(data[[1]]$coefs$A_d)[3]
  
  c_tvp <- unlist(lapply(data, function(x){return(x$specs$tvp[1])}))
  tvp <- any(c_tvp)
  c_structural <- unlist(lapply(data, function(x){return(x$specs$structural)}))
  structural <- any(c_structural)
  c_tvp_struct <- unlist(lapply(data, function(x){return(x$specs$tvp[2])}))
  tvp_struct <- any(c_tvp_struct)
  c_sv <- unlist(lapply(data, function(x){return(x$specs$sv)}))
  sv <- any(c_sv)
  t.total <- max(unlist(lapply(data, function(x){return(dim(x$coefs$A_d)[2])})))
  
  const.w <- any(unlist(lapply(W, function(x){is.na(dim(x)[3])})))
  t.final <- 1
  if (is.null(t)){
    if (tvp | !const.w) {
      t.final <- max(t.total, unlist(lapply(W, function(x){dim(x)[3]})) - data[[1]]$specs$lags$maximum.lag)
    }
    if (const.w) {
      t.final <- t.total 
    }
  }
  
  t.sv <- 1
  if (sv & is.null(t)) {
    t.sv <- max(unlist(lapply(data, function(x){return(dim(x$coefs$Omega)[2])}))) 
  }
  Omega <- array(0, c(k * k, t.sv, draws))

  if (structural) {
    if (tvp_struct) {
      A0 <- array(NA, c(k * k, t.final, draws)) 
    } else {
      A0 <- array(NA, c(k * k, 1, draws))
    }
    if (sv) {
      Sigma <- array(NA, c(k * k, t.sv, draws)) 
    } else {
      Sigma <- array(NA, c(k * k, 1, draws)) 
    }
  }
  
  G0_i <- array(NA, c(k * k, t.final, draws))
  
  if (p > 0){
    G <- array(NA, c(k * k * p, t.final, draws))
  } else {
    G <- NULL
  }
  
  if (global){
    n.global <- max(n.g)
    G_g <- array(0, c(k * n.global * (p.global + 1), t.final, draws))
  } else {
    G_g <- NULL
  }
  
  det <- unique(unlist(lapply(data, function(x) {x$specs$deterministic.terms}))) != "I"
  if (det){
    n.det <- unlist(lapply(data, function(x){return(dim(x$coefs$A_det)[1] / dim(x$x)[2])}))
    G_det <- array(0, c(k * max(n.det), t.final, draws))
  } else {
    G_det <- NULL
  }
  
  Omega_temp <- matrix(0, k, k)
  
  if (is.null(t) & (tvp | !const.w)) {
    t.final <- 1:t.final
  }
  if (is.null(t) & sv) {
    t.sv <- 1:t.sv
  }
  
  t.par <- 1
  t.choice <- NULL
  if (!is.null(t)) {
    t.choice <- t
    if (tvp) {
      t.par <- t
    }
    if (!const.w) {
      t.w <- t
    }
  }
  
  cat("Solving GVAR model.\n")
  pb <- utils::txtProgressBar(width = 70, style = 3)
  for (draw in 1:draws){
    for (t in t.final){
      if (is.null(t.choice)) {
        if (tvp){
          t.par <- t
        }
        if (!const.w) {
          t.w <- t
        }
      }
      
      G0_temp <- matrix(NA, k, k)
      G_temp <- matrix(0, k, k * p)
      if (global){
        G_g_temp <- matrix(NA, k, n.global * (p.global + 1))
      }
      if (det){
        G_det_temp <- matrix(NA, k, max(n.det))
      }
      
      # Contemporaneous effects
      for (i in 1:n.c){
        t_d <- ifelse(dim(data[[i]]$coefs$A_d)[2] == 1, t, t.par)
        pos.ki <- pos.i[i] + 1:k.i[i]
        lag.i <- max(c(p.i[i], p.star.i[i]))
        
        if (const.w) {
          G0_temp[pos.ki, ] <- cbind(diag(1, k.i[i]),
                                     -matrix(data[[i]]$coefs$A_s_0[, t.par, draw], k.i[i]))%*%W[[i]]
        } else {
          G0_temp[pos.ki, ] <- cbind(diag(1, k.i[i]),
                                     -matrix(data[[i]]$coefs$A_s_0[, t.par, draw], k.i[i]))%*%W[[i]][,, t.w + data[[i]]$specs$lags$maximum.lag]
        }
        
        # Lags
        if (lag.i > 0){
          for (j in 1:lag.i){
            if (j <= p.i[i]){
              temp.d <- matrix(data[[i]]$coefs$A_d[(j - 1) * (k.i[i]^2) + 1:(k.i[i]^2), t_d, draw], k.i[i])
            } else {
              temp.d <- matrix(0, k.i[i], k.i[i])  
            }
            
            if (j <= p.star.i[i]){
              temp.s <- matrix(data[[i]]$coefs$A_s[(j - 1) * k.i[i] * k.s.i[i] + 1:(k.i[i] * k.s.i[i]), t.par, draw], k.i[i])
            } else {
              temp.s <- matrix(0, k.i[i], k.s.i[i])
            }
            
            if (const.w){
              G_temp[pos.ki, (j - 1) * k + 1:k] <- cbind(temp.d, temp.s) %*% W[[i]]
            } else {
              G_temp[pos.ki, (j - 1) * k + 1:k] <- cbind(temp.d, temp.s) %*% W[[i]][,, t.w + data[[i]]$specs$lags$maximum.lag]
            }
            rm(list = c("temp.d", "temp.s"))
          }
        }
        
        # Global variables
        if (global & !is.na(p.global.i[i])){
          if (p.global.i[i] > 0){
            for (j in 1:(p.global.i[i] + 1)){
              G_g_temp[pos.ki, (j - 1) * n.global + 1:n.global] <- data[[i]]$coefs$A_g[(j - 1) * k.i[i] * n.g[i] + 1:(k.i[i] * n.g[i]), t.par, draw]
            }
          }   
        }
        
        # Deterministic terms
        if (det){
          G_det_temp[pos.ki,] <- matrix(data[[i]]$coefs$A_det[, t.par, draw], k.i[i])
        }
      }
      
      G0_i_temp <- solve(G0_temp)
      G0_i[, t, draw] <- G0_i_temp
      G[, t, draw] <- G0_i_temp %*% G_temp
      if (global){
        G_g[, t, draw] <- G0_i_temp %*% G_g_temp
      }
      if (det){
        G_det[, t, draw] <- G0_i_temp %*% G_det_temp
      }
    }
    
    if (structural) {
      A0_temp <- matrix(0, k, k)
      for (i in 1:n.c){
        pos.ki <- pos.i[i] + 1:k.i[i]
        A0_temp[pos.ki, pos.ki] <- matrix(data[[i]]$coefs$A0[, 1, draw], k.i[i]) 
      }
      A0[, t, draw] <- A0_temp
      rm(A0_temp)
      
      Sigma_temp <- matrix(0, k, k)
      for (i in 1:n.c){
        pos.ki <- pos.i[i] + 1:k.i[i]
        diag(Sigma_temp)[pos.ki] <- diag(matrix(data[[i]]$coefs$Sigma[, 1, draw], k.i[i]))
      }
      Sigma[, t, draw] <- Sigma_temp
      rm(Sigma_temp)
    }
    
    if (sv){
      for (t in t.sv){
        if (is.null(t.choice)) {
          t.s <- t
        } else {
          t.s <- t.choice
        }
        for (i in 1:n.c){
          pos.ki <- pos.i[i] + 1:k.i[i]
          Omega_temp[pos.ki, pos.ki] <- matrix(data[[i]]$coefs$Omega[, t.s, draw], k.i[i]) 
        }
        Omega[, t, draw] <- Omega_temp
      }
    } else {
      for (i in 1:n.c){
        pos.ki <- pos.i[i] + 1:k.i[i]
        Omega_temp[pos.ki, pos.ki] <- matrix(data[[i]]$coefs$Omega[, 1, draw], k.i[i]) 
      }
      Omega[, 1, draw] <- Omega_temp
    }
    
    utils::setTxtProgressBar(pb, draw/draws)
  }
  close(pb)
  
  gvar_coefs <- list(G0_i = G0_i, G = G, G_det = G_det, G_g = G_g, Omega = Omega)
  
  if (structural) {
    gvar_coefs$A0 <- A0
    gvar_coefs$Sigma <- Sigma
  }
  
  structural <- unlist(lapply(data, function(x) {return(x$specs$structural)}))
  
  result <- list(country.models = original.data,
                 GVAR = list(data = list(X = global.data$X, X.global = global.data$X.global),
                             coefs = gvar_coefs,
                             specs = list(lags = c("endogenous" = p, "global" = p.global),
                                          type = data[[1]]$specs$type,
                                          case = data[[1]]$specs$deterministic.terms,
                                          structural = structural,
                                          tvp = c_tvp,
                                          sv = c_sv,
                                          index = index,
                                          criteria = criteria)))
  return(result)
}