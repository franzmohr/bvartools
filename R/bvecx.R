#' Estimate a Vector Error Correction Model
#' 
#' Estimates a country-specific VECX* model with constant parameters.
#' 
#' @param data A list containing the data necessary to estimate a GVAR country model.
#' @param iterations Number of iterations of the Gibbs sampler.
#' @param burnin Number of burnin draws of the Gibbs sampler.
#' @param thin Thinning factor.
#' 
#' @return Some object.
#' 
#' @export
bvecx <- function(data, iterations = 15000, burnin = 5000, thin = NULL){
  estimation.data <- gen.vecx(data)
  y <- estimation.data$y
  ect <- estimation.data$ect
  x <- estimation.data$x
  t <- dim(y)[2]
  global <- !is.null(data$x.g)
  
  r <- estimation.data$r
  if (is.na(r)) {stop("Rank of Pi must be specified.")}
  if (r > 0) {ecm <- TRUE} else {ecm <- FALSE}
  
  # Estimate the model
  store <- iterations-burnin
  if (store > (iterations - burnin)) {stop("Number of iterations must be larger than burn-in draws.")}
  
  n <- dim(y)[1]
  n.ect <- dim(ect)[1]
  n.x <- dim(x)[1]
  n.B <- n*n.x
  n.det.r <- estimation.data$deterministic["restricted"]
  n.det.ur <- estimation.data$deterministic["unrestricted"]
  
  n.alpha <- r*n
  n.alpha.B <- n.alpha + n.B
  n.beta <- r*n.ect
  
  prior <- data$priors
  B.mu.prior <- prior$B$constant[[1]]
  B.V.i.prior <- prior$B$constant[[2]]
  Sigma.df.post <- t#+prior$Sigma$constant[[1]]
  Sigma.V.prior <- prior$Sigma$constant[[2]]
  
  shrinkage <- prior$Shrinkage$type != "none"
  if (shrinkage){
    shrinkage.type <- prior$Shrinkage$type
    if (shrinkage.type == "BVS"){
      Gamma.B <- matrix(1, n.alpha.B)
      restricted.variables <- prior$Shrinkage$spec$B[[1]]
      B.lpr.include.prior <- prior$Shrinkage$spec$B[[2]]
      B.lpr.exclude.prior <- prior$Shrinkage$spec$B[[3]]
    }
  }
  
  x.B <- rbind(matrix(NA, r, t), x)
  beta <- matrix(stats::rnorm(n.beta), n.ect, r)
  
  Sigma <- diag(.1, n)
  Sigma.i <- solve(Sigma)
  
  y.beta <- y
  
  draws.B <- matrix(NA, n.B, iterations - burnin)
  draws.Pi <- matrix(NA, n*n.ect, iterations - burnin)
  if (r == 0) {draws.Pi[,] <- 0}
  draws.Sigma <- matrix(NA, n*n, iterations - burnin)
  draws.LL <- matrix(NA, t, iterations - burnin)
  
  if (iterations <= 100){
    update <- 1:iterations
  } else {
    update <- round(seq(from = 0, to = iterations, length.out = 100)) # Refreshment rate of the progress bar
  }
  pb <- utils::txtProgressBar(width = 70, style = 3)
  
  for (draw in 1:iterations){
    if (ecm){
      x.B[1:r,] <- crossprod(beta, ect)
    }
    if (shrinkage) {
      if (shrinkage.type == "BVS"){
        Z.B <- kronecker(t(x.B), diag(1, n))
        Z.B.res <- Z.B%*%diag(Gamma.B[,1], n.alpha.B)
        B <- sur_normal(y, Z.B.res, Sigma.i, B.mu.prior, B.V.i.prior)
        Gamma.B <- bvs(y, Z.B, B, Gamma.B, Sigma.i, restricted.variables, B.lpr.include.prior, B.lpr.exclude.prior)
        B <- diag(Gamma.B[, 1], n.alpha.B)%*%B
      }
    } else {
      B <- minn_normal(y, x.B, Sigma.i, B.mu.prior, B.V.i.prior)
    }
    B <- matrix(B, n)
    if (ecm){
      alpha <- matrix(B[,1:r], n)
      if (n.x > 0){
        B <- matrix(B[, -(1:r)], n)
        y.beta <- y - B%*%x 
      }
      Pi.temp <- minn_Pi(y.beta, ect, alpha, Sigma.i)
      Pi <- Pi.temp[[1]]
      beta <- Pi.temp[[2]]
      y.S <- y.beta - Pi%*%ect
    } else {
      y.S <- y - B%*%x
    }
    
    # Error covariance matrix
    Sigma.i <- stats::rWishart(1, Sigma.df.post, solve(tcrossprod(y.S)))[,,1]
    Sigma <- solve(Sigma.i)
    
    # Likelihood calculations
    if (draw>burnin){
      pos.draw <- draw - burnin
      if (n.x>0){
        draws.B[, pos.draw] <- B
      }
      if (ecm){
        draws.Pi[, pos.draw] <- Pi
      }
      draws.Sigma[, pos.draw] <- Sigma
      draws.LL[, pos.draw] <- getLL(y.S, Sigma, Sigma.i)
    }
    # Update progress bar
    if (is.element(draw,update)) {utils::setTxtProgressBar(pb, draw/iterations)}
  }
  close(pb)
  
  # Transform VECX parameters to VARX
  names.y <- gsub("d.", "", dimnames(y)[[1]])
  names.ect <- dimnames(ect)[[1]]
  names.x <- gsub("d.", "", dimnames(x)[[1]])
  
  if (is.null(thin)){
    tw <- seq(to = store, length.out = store)
  } else {
    tw <- sample(x = 1:store, size = store/thin, replace = FALSE)
  }
  
  draws <- length(tw)
  A0 <- array(diag(1, n), dim = c(n*n, 1, draws))
  dimnames(A0) <- list(paste(names.y, rep(names.y, each = n), sep = "_"), NULL, NULL)
  B <- array(draws.B[, tw], dim = c(n.B, 1, draws))
  Pi <- array(draws.Pi[, tw], dim = c(n*n.ect, 1, draws))
  Sigma <- array(draws.Sigma[, tw], dim = c(n*n, 1, draws))
  dimnames(Sigma) <- list(paste(names.y, rep(names.y, each = n), sep = "_"), NULL, NULL)
  
  # Domestic variables
  p <- estimation.data$domestic["lag"]
  n.d <- estimation.data$domestic["total"]
  if (p > 0){
    B.d <- array(B[1:(n*n.d),,], c(n*n.d, 1, draws))
  }
  
  p <- p + 1
  if (p > 1){
    W <- diag(-1, n*p)
    W[1:n, 1:n] <- diag(1, n)
    W[-(1:n), -(n*(p-1) + 1:n)] <- W[-(1:n), -(n*(p - 1) + 1:n)] + diag(n*(p - 1))
    J <- matrix(0, n, n*p)
    J[1:n, 1:n] <- diag(1, n)
    
    A.d <- array(NA,c(n*n*p, 1, draws))
    for (draw in 1:draws){
      A.d[, 1, draw] <- cbind(matrix(Pi[1:(n*n), 1, draw], n), matrix(B.d[, 1, draw], n))%*%W + J
    }
    
    names.d <- NULL
    for (i in 1:p){
      names.d <- c(names.d, paste(names.y, ".l", i, sep = ""))
    }
    dimnames(A.d) <- list(paste(names.y, rep(names.d, each = n), sep = "_"), NULL, NULL)
  } else {
    A.d <- array(NA, c(n*n, 1, draws))
    for (draw in 1:draws){
      A.d[, 1, draw] <- matrix(Pi[1:(n*n), 1, draw], n) + diag(1, n)
    }
    dimnames(A.d) <- list(paste(names.y, rep(names.y, each = n), sep = "_"), NULL, NULL)
  }
  
  # Foreign variables
  n.star <- estimation.data$foreign["dim"]
  q <- estimation.data$foreign["lag"]
  n.s <- estimation.data$foreign["total"]
  B.s <- array(B[n.d*n + 1:(n*n.s),,], c(n*n.s, 1, draws))
  
  q <- q + 1
  if (q > 0){
    W.s <- diag(-1, n.star*(q + 1))
    W.s[1:n.star, 1:n.star] <- 0
    W.s[1:n.star, n.star + 1:n.star] <- diag(1, n.star)
    W.s[-(1:n.star), 1:(n.star*q)] <- W.s[-(1:n.star), 1:(n.star*q)] + diag(1, n.star*q)
    
    A.s <- array(NA, c(n*n.star*(q + 1), 1, draws))
    for (draw in 1:draws){
      A.s[, 1, draw] <- cbind(matrix(Pi[n*n + 1:(n*n.star), 1, draw], n), matrix(B.s[, 1, draw], n))%*%W.s
    }
    
    names.s <- names.x[n.d + 1:n.star]
    for (i in 1:q){
      names.s <- c(names.s, paste(names.x[n.d + 1:n.star], ".l", i, sep = ""))
    }
    dimnames(A.s) <- list(paste(names.y, rep(names.s, each = n), sep = "_"), NULL, NULL)
  } else {
    stop("Lag of foreign star variables in VARX model must be at least 1.")
  }
  
  A.s.0 <- array(A.s[1:(n.star*n), 1,], c(n*n.star, 1, draws))
  dimnames(A.s.0) <- list(dimnames(A.s)[[1]][1:(n*n.star)], NULL, NULL)
  names.A.s <- dimnames(A.s)[[1]][-(1:(n*n.star))]
  A.s <- array(A.s[-(1:(n*n.star)), 1,], c(n*n.s, 1, draws))
  dimnames(A.s)[[1]] <- names.A.s
  
  # Global variables
  A.g <- NULL
  n.global <- estimation.data$global["dim"]
  s <- estimation.data$global["lag"]
  n.g <- estimation.data$global["total"]
  if (global){
    if (!is.na(s)){
      names.g <- gsub("d.", "", dimnames(x)[[1]][n.d + n.s + 1:n.global])
      B.g <- array(B[n*n.d + n*n.s + 1:(n*n.g), 1,], c(n*n.g, 1, draws))
      A.g <- array(NA,c(n*(n.g + 1), 1, draws))
      if (s > 1){
        W.g <- diag(-1, n.global*(s + 1))
        W.g[1:n.global, 1:n.global] <- 0
        W.g[1:n.global, n.global + 1:n.global] <- diag(1, n.global)
        W.g[-(1:n.global), 1:(n.global*s)] <- W.g[-(1:n.global), 1:(n.global*s)] + diag(1, n.global*s)
        
        for (draw in 1:draws){
          A.g[, 1, draw] <- cbind(matrix(Pi[n*(n + n.star) + 1:(n*n.global), 1, draw], n), matrix(B.g[, 1, draw], n))%*%W.g
        }
        
        names.global <- names.g
        for (i in 1:s) {
          names.global <- c(names.global, paste(names.g, ".l", i, sep = ""))
        }
        dimnames(A.g) <- list(paste(names.y, rep(names.global, each = n), sep="_"), NULL, NULL)
      } else {
        for (draw in 1:draws){
          A.g[, 1, draw] <- cbind(matrix(B.g[, 1, draw], n),
                                  matrix(Pi[n*(n + n.star) + 1:(n*n.global), 1, draw], n) - matrix(B.g[, 1, draw], n))
        }
        names.global <- c(names.x[n.d + n.s + 1:n.global], paste(names.x[n.d + n.s + 1:n.global], ".l", 1, sep = ""))
        dimnames(A.g) <- list(paste(names.y, rep(names.global, each = n), sep = "_"), NULL, NULL)
      }
    }
  }
  
  n.det <- n.det.r + n.det.ur
  det <- array(NA, c(n.det*n, 1, draws))
  names.det <- NULL
  # Unrestricted determinisitc terms
  if (n.det.ur > 0){
    det[1:(n*n.det.ur), 1,] <- B[n*(n.d + n.s + n.g) + 1:(n.det.ur*n), 1,]
    names.det <- c(names.det, names.x[n.d + n.s + n.g + 1:n.det.ur])
  }
  # Restriced deterministic terms
  if (n.det.r > 0){
    det[(n*n.det.ur) + 1:(n*n.det.r), 1,] <- Pi[n*(n + n.star + n.global) + 1:(n*n.det.r), 1,]
    names.det <- c(names.det, names.ect[n + n.star + n.global + 1:n.det.r])
  }
  if (!is.null(det)){
    dimnames(det) <- list(paste(names.y, rep(names.det, each = n), sep = "_"), NULL, NULL)
  }
  
  # Test statistic
  test.stat <- c("AIC" = 2*sum(log(rowMeans(draws.LL))) - (n.alpha.B)*2,
                 "BIC" = 2*sum(log(rowMeans(draws.LL))) - (n.alpha.B)*log(t),
                 "HQ" = 2*sum(log(rowMeans(draws.LL))) - (n.alpha.B)*2*log(log(t)))
  
  # Prepare output object
  result <- c(data,list(coefs = list(A0 = A0, A.d = A.d, A.s.0 = A.s.0, A.s = A.s, A.g = A.g, Det = det, Sigma = Sigma),
                        criteria = list("Criteria" = test.stat)))
  return(result)
}