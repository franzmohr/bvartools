#' Estimate a Vector Autoregressive Model
#' 
#' Estimates a country-specific VARX* model with constant parameters.
#' 
#' @param data A list containing the data necessary to estimate a GVAR country model.
#' @param iterations Number of iterations of the Gibbs sampler.
#' @param burnin Number of burnin draws of the Gibbs sampler.
#' @param thin Thinning factor.
#' 
#' @return Some object.
#' 
#' @export
bvarx <- function(data, iterations = 15000, burnin = 5000, thin = NULL){
  estimation.data <- gen.varx(data)
  y <- estimation.data$y
  x <- estimation.data$x
  t <- dim(y)[2]
  global <- !is.null(data$x.g)
  
  store <- iterations - burnin
  if (store > (iterations - burnin)){stop("Number of iterations must be larger than burn-in draws.")}
  
  n <- dim(y)[1]
  n.x <- dim(x)[1]
  n.A <- n * n.x
  n.det <- estimation.data$deterministic
  
  prior <- data$priors
  A.mu.prior <- prior$A$constant[[1]]
  A.V.i.prior <- prior$A$constant[[2]]
  Sigma.df.post <- t + prior$Sigma$constant[[1]]
  Sigma.V.prior <- prior$Sigma$constant[[2]]
  
  shrinkage <- prior$Shrinkage$type != "none"
  if (shrinkage){
    shrinkage.type <- prior$Shrinkage$type
    if (shrinkage.type == "BVS"){
      Gamma.A <- matrix(1, n.A)
      Z.A <- kronecker(t(x), diag(1, n))
      restricted.variables <- prior$Shrinkage$spec$A[[1]]
      A.lpr.include.prior <- prior$Shrinkage$spec$A[[2]]
      A.lpr.exclude.prior <- prior$Shrinkage$spec$A[[3]]
    }
  }
  
  Sigma <- diag(.0001, n)
  Sigma.i <- solve(Sigma)
  
  draws.A <- matrix(NA, n.A, iterations - burnin)
  draws.Sigma <- matrix(NA, n*n, iterations - burnin)
  draws.LL <- matrix(NA, t, iterations - burnin)
  
  if (iterations <= 100){
    update <- 1:iterations
  } else {
    update <- round(seq(from = 0, to = iterations, length.out = 100)) # Refreshment rate of the progress bar
  }
  pb <- utils::txtProgressBar(width = 70, style = 3) # Create progress bar
  
  for (draw in 1:iterations){
    if (shrinkage){
      if (shrinkage.type == "BVS"){
        Z.A.res <- Z.A%*%diag(Gamma.A[, 1], n.A)
        A <- posterior_normal_sur(y, Z.A.res, Sigma.i, A.mu.prior, A.V.i.prior)
        
        bla <- matrix(0, n.A, n.A)
        for (i in 1:t) {
          bla <- bla + t(Z.A.res[(i-1)*n + 1:n,]) %*% Sigma.i %*% Z.A.res[(i-1)*n + 1:n,]
        }
        
        
        Gamma.A <- bvs(y, Z.A, A, Gamma.A, Sigma.i, restricted.variables, A.lpr.include.prior, A.lpr.exclude.prior)
        A <- diag(Gamma.A[, 1], n.A)%*%A
      }
    } else {
      A <- posterior_normal(y, x, Sigma.i, A.mu.prior, A.V.i.prior)
    }
    
    y.temp <- y - matrix(A, n)%*%x
    
    # Error covariance matrix
    Sigma.i <- stats::rWishart(1, Sigma.df.post, solve(Sigma.V.prior + tcrossprod(y.temp)))[,,1]
    Sigma <- solve(Sigma.i)
    
    # Likelihood calculations
    if (draw > burnin){
      pos.draw <- draw - burnin
      draws.A[, pos.draw] <- A
      draws.Sigma[, pos.draw] <- Sigma
      draws.LL[, pos.draw] <- getLL(y.temp, Sigma, Sigma.i)
    }
    # Update progress bar
    if (is.element(draw, update)) {utils::setTxtProgressBar(pb, draw/iterations)}
  }
  close(pb)
  
  names.y <- dimnames(y)[[1]]
  names.x <- dimnames(x)[[1]]
  
  p <- estimation.data$domestic["lag"]
  q <- estimation.data$foreign["lag"]
  n.star <- estimation.data$foreign["dim"]
  n.s <- estimation.data$foreign["total"]
  n.g <- estimation.data$global["total"]
  
  if (is.null(thin)){
    tw <- seq(to = store, length.out = store)
  } else {
    tw <- sample(x = 1:store, size = store/thin, replace = FALSE)
  }
  
  draws <- length(tw)
  A0 <- array(diag(1, n), dim = c(n*n, 1, draws))
  dimnames(A0) <- list(paste(names.y, rep(names.y, each = n), sep = "_"), NULL, NULL)
  A.d <- array(draws.A[1:(n*n*p), tw], dim = c(n*n*p, 1, draws))
  dimnames(A.d) <- list(paste(names.y, rep(names.x[1:(n*p)], each = n), sep = "_"), NULL, NULL)
  A.s.0 <- array(draws.A[n*n*p+1:(n*n.star), tw], dim = c(n*n.star, 1, draws))
  dimnames(A.s.0) <- list(paste(names.y, rep(names.x[n*p+1:n.star], each = n), sep = "_"), NULL, NULL)
  A.s <- array(draws.A[n*(n*p + n.star) + 1:(q*n*n.star), tw], dim = c(q*n*n.star, 1, draws))
  dimnames(A.s) <- list(paste(names.y, rep(names.x[n*p + n.star + 1:(n.star*q)], each = n), sep = "_"), NULL, NULL)
  if (global){
    s <- estimation.data$global["lag"]
    A.g <- array(draws.A[n*(n*p + n.star*(1 + q)) + 1:(n*n.g), tw], dim = c(n*n.g, 1, draws))
    dimnames(A.g) <- list(paste(names.y, rep(names.x[n*p + n.s + 1:n.g], each = n), sep = "_"), NULL, NULL)
  } else {
    A.g <- NULL
  }
  if (n.det > 0){
    det <- array(draws.A[n*(n*p + n.star*(1 + q) + n.g) + 1:(n*n.det), tw], dim = c(n*n.det, 1, draws))
    dimnames(det) <- list(paste(names.y, rep(names.x[n*p + n.s + n.g + 1:n.det], each = n), sep = "_"), NULL, NULL)
  } else {
    det <- NULL
  }
  Sigma <- array(draws.Sigma[, tw], dim = c(n*n, 1, draws))
  dimnames(Sigma) <- list(paste(names.y, rep(names.y, each = n), sep = "_"), NULL, NULL)
  
  test.stat <- c("AIC" = 2*sum(log(rowMeans(draws.LL))) - n.A*2,
                 "BIC" = 2*sum(log(rowMeans(draws.LL))) - n.A*log(t),
                 "HQ" = 2*sum(log(rowMeans(draws.LL))) - n.A*2*log(log(t)))
  
  result <- c(data, list(coefs = list(A0 = A0, A.d = A.d, A.s.0 = A.s.0, A.s = A.s, A.g = A.g, Det = det, Sigma = Sigma),
                         criteria = list("Criteria" = test.stat)))
  return(result) 
}
