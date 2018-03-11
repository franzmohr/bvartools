#' Solve a GVAR model
#' 
#' Combines the country model results to a global VAR model and solves it.
#' 
#' @param data a list containing the results of the country estimates.
#' @param global.data  a list as produced by the function \code{\link{country_models}}.
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
solve_gvar <- function(data, global.data, ic = "BIC", t = NULL){
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
  if (ic%in%c("AIC", "BIC", "HQ")){
    for (i in names(data)){
      c.pos <- which(criteria[, "Country"] == i)
      pos.temp <- which.max(criteria[c.pos, ic])
      criteria[c.pos[pos.temp], "Maximum"] <- TRUE
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
  tvp <- unique(unlist(lapply(data, function(x){return(dim(x$coefs$A0)[2])}))) > 1
  t.sv <- unique(unlist(lapply(data, function(x){return(dim(x$coefs$Sigma)[2])})))
  sv <- all(t.sv > 1)
  const.w <- any(unlist(lapply(W, function(x){is.na(dim(x)[3])})))
  
  # Determine max numbers of lags
  n.c <- length(data)
  k.i <- unlist(lapply(data, function(x){length(x$specs$domestic.variables)}))
  k.s.i <- unlist(lapply(data, function(x){length(x$specs$foreign.variables)}))
  p.i <- unlist(lapply(data, function(x){x$specs$lags$lags$domestic}))
  q.i <- unlist(lapply(data, function(x){x$specs$lags$lags$foreign}))
  n.g <- unlist(lapply(data, function(x){g <- x$specs$global.variables
  if (is.na(g)){return(0)} else {return(length(g))}
  }))
  s.i <- unlist(lapply(data, function(x){x$specs$lags$lags$global}))
  
  p <- max(p.i, q.i)
  if (all(is.na(s.i))){
    s <- NA
  } else {
    s <- max(stats::na.omit(s.i))
  }
  global <- !is.na(s)
  
  pos.i <- c(0, cumsum(k.i))
  
  draws <- dim(data[[1]]$coefs$A0)[3]

    t.total <- max(unlist(lapply(data, function(x){return(dim(x$coefs$A0)[2])})))
    if (!const.w) {
      t.final <- max(t.total, unlist(lapply(W, function(x){dim(x)[3]})))
    } else {
      t.final <- t.total
    } 
  
  G0.i <- array(NA, c(k * k, t.final, draws))
  if (p > 0){
    G <- array(NA, c(k * k * p, t.final, draws))
  } else {
    G <- NULL
  }
  
  if (global){
    n.global <- max(n.g)
    G.g <- array(0, c(k * n.global * (s + 1), t.final, draws))
  } else {
    G.g <- NULL
  }
  
  det <- unique(unlist(lapply(data, function(x) {x$specs$deterministic.terms}))) != "I"
  if (det){
    n.det <- unlist(lapply(data, function(x){return(dim(x$coefs$Det)[1]/dim(x$x)[2])}))
    G.det <- array(0, c(k * max(n.det), t.final, draws))
  } else {
    G.det <- NULL
  }
  
  Sigma <- array(0, c(k * k, t.sv, draws))
  Sigma.temp <- matrix(0, k, k)
  
  cat("Solving GVAR model.\n")
  pb <- utils::txtProgressBar(width = 70, style = 3)
  for (draw in 1:draws){
    for (tt in 1:t.final){
      if (const.w) {
        t <- tt
      } else {
        if (tvp) {
          t <- tt
        } else {
          t <- 1
        }
      }
      G0.temp <- matrix(NA, k, k)
      G.temp <- matrix(0, k, k * p)
      if (global){
        G.g.temp <- matrix(NA, k, n.global * (s + 1))
      }
      if (det){
        G.det.temp <- matrix(NA, k, max(n.det))
      }
      
      for (i in 1:n.c){
        pos.ki <- pos.i[i] + 1:k.i[i]
        lag.i <- max(c(p.i[i], q.i[i]))
        
        if (const.w) {
          G0.temp[pos.ki, ] <- cbind(matrix(data[[i]]$coefs$A0[, t, draw], k.i[i]),
                                     -matrix(data[[i]]$coefs$A.s.0[, t, draw], k.i[i]))%*%W[[i]]
        } else {
          G0.temp[pos.ki, ] <- cbind(matrix(data[[i]]$coefs$A0[, t, draw], k.i[i]),
                                     -matrix(data[[i]]$coefs$A.s.0[, t, draw], k.i[i]))%*%W[[i]][,,tt]
        }
        
        if (lag.i > 0){
          for (j in 1:lag.i){
            if (j <= p.i[i]){
              temp.d <- matrix(data[[i]]$coefs$A.d[(j - 1)*(k.i[i]^2) + 1:(k.i[i]^2), t, draw], k.i[i])
            } else {
              temp.d <- matrix(0, k.i[i], k.i[i])  
            }
            
            if (j <= q.i[i]){
              temp.s <- matrix(data[[i]]$coefs$A.s[(j - 1)*k.i[i]*k.s.i[i] + 1:(k.i[i]*k.s.i[i]), t, draw], k.i[i])
            } else {
              temp.s <- matrix(0, k.i[i], k.s.i[i])
            }
            
            if (const.w){
              G.temp[pos.ki, (j-1)*k + 1:k] <- cbind(temp.d, temp.s)%*%W[[i]]
            } else {
              G.temp[pos.ki, (j-1)*k + 1:k] <- cbind(temp.d, temp.s)%*%W[[i]][,,tt]
            }
            rm(list = c("temp.d", "temp.s"))
          }
        }
        
        if (global & !is.na(s.i[i])){
          if (s.i[i] > 0){
            for (j in 1:(s.i[i]+1)){
              G.g.temp[pos.ki, (j-1)*n.global + 1:n.global] <- data[[i]]$coefs$A.g[(j - 1)*k.i[i]*n.g[i] + 1:(k.i[i]*n.g[i]), t, draw]
            }
          }   
        }
        
        if (det){
          G.det.temp[pos.ki,] <- matrix(data[[i]]$coefs$Det[, t, draw], k.i[i])
        }
      }
      
      G0.i.temp <- solve(G0.temp)
      G0.i[, tt, draw] <- G0.i.temp
      G[, tt, draw] <- G0.i.temp%*%G.temp
      if (global){
        G.g[, tt, draw] <- G0.i.temp%*%G.g.temp
      }
      if (det){
        G.det[, tt, draw] <- G0.i.temp%*%G.det.temp
      }
    }
    if (sv){
      for (t in 1:t.sv){
        for (i in 1:n.c){
          pos.ki <- pos.i[i]+1:k.i[i]
          Sigma.temp[pos.ki, pos.ki] <- matrix(data[[i]]$coefs$Sigma[, t, draw], k.i[i]) 
        }
        Sigma[, t, draw] <- Sigma.temp
      }
    } else {
      for (i in 1:n.c){
        pos.ki <- pos.i[i] + 1:k.i[i]
        Sigma.temp[pos.ki, pos.ki] <- matrix(data[[i]]$coefs$Sigma[, 1, draw], k.i[i]) 
      }
      Sigma[,1,draw] <- Sigma.temp
    }
    utils::setTxtProgressBar(pb, draw/draws)
  }
  close(pb)
  
  gvar.coefs <- list(G0.i = G0.i, G = G, G.det = G.det, G.g = G.g, Sigma = Sigma,
                     lags = c("endogenous" = p, "global" = s))
  
  structural <- unlist(lapply(data, function(x) {return(x$specs$structural)}))
  tvp <- unlist(lapply(data, function(x) {return(x$specs$tvp)}))
  sv <- unlist(lapply(data, function(x) {return(x$specs$sv)}))
  
  result <- list(Country.models = original.data,
                 GVAR = list(data = list(X = global.data$X, X.global = global.data$X.global),
                             coefs = gvar.coefs,
                             specs = list(type = type,
                                          case = data[[1]]$specs$deterministic.terms,
                                          structural = structural,
                                          tvp = tvp,
                                          sv = sv,
                                          index = index,
                                          criteria = criteria)))
  return(result)
}