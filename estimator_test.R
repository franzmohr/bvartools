rm(list=ls())

################### Load functions ###################
devtools::load_all(".")

set.seed(1)

# Generate the VEC model
tt <- 1000
n <- 3
n.x <- 2
n.ect <- n + n.x

beta <- matrix(c(0, 1, 0, 2, 0,
                 0, 0, 1, 0, -4), n.ect * 2, tt)
for (i in 2:tt) {
  beta[10, i] <- beta[10, i] + i * .01
}
alpha <- matrix(c(0, -.1, 0, 0, 0, -.05), n)
Pi <- matrix(NA, n * n.ect, tt)
for (i in 1:tt) {
  Pi[, i] <- alpha %*% t(matrix(beta[, i], n.ect))
}

A0 <- diag(1, n)
A0[lower.tri(A0)] <- c(.03, .4, -.1)
A0_i <- solve(A0)
Sigma <- diag(sqrt(.1), n)

A1 <- diag(.3, n)
B0 <- diag(.5, n)[, -1]
B1 <- diag(.2, n)[, -3]

C1 <- array(NA, c(n, n, tt))
for (i in 1:tt) {
  C1[,,i] <- matrix(Pi[,i], n)[, 1:n] + diag(1, n) + A1 
}
C2 <- -A1
D0 <- B0
D1 <- array(NA, c(n, 2, tt))
for (i in 1:tt) {
  D1[,,i] <- matrix(Pi[,i], n)[, -(1:n)] + B1 - B0 
}
D2 <- -B1

Pi_true <- Pi
#true <- cbind(Pi, A1, B0, B1)

#x <- diag(c(.01, .02), 2)  %*% t(matrix(1:(tt + 2), tt + 2, 1))# + matrix(rnorm(2 * (tt + 2)), 2)
x <- t(apply(matrix(rnorm(2 * (tt + 2)), 2), 1, cumsum))
y <- matrix(0, n, tt + 2)

for (i in 2 + 1:tt) {
  y[, i] <- C1[,,i - 2] %*% y[, i - 1] + C2 %*% y[, i - 2] + D0 %*% x[, i] + D1[,,i - 2] %*% x[, i - 1] + D2 %*% x[, i - 2] + A0_i %*% Sigma %*% rnorm(n)
}

pos <- 2:203
y <- ts(t(y[, pos])); dimnames(y)[[2]] <- c("y1", "y2", "y3")
x <- ts(t(x[, pos])); dimnames(x)[[2]] <- c("x1", "x2")
plot.ts(cbind(y, x))

#####################################################################################
# Generate data for VEC
level.domestic <- y
level.star <- x
ect <- cbind(level.domestic, level.star)
n.ect <- dim(ect)[2]
ect <- stats::lag(ect, -1)

diff.domestic <- diff(level.domestic)
total <- cbind(diff.domestic, ect)
# Lags of domestic variables
for (i in 1:1){
  total <- cbind(total, stats::lag(diff.domestic, -i))
}

# Lags of star variables
diff.star <- diff(x)
total <- cbind(total, diff.star)
for (i in 1:1){
  total <- cbind(total, stats::lag(diff.star, -i))
}

total <- stats::na.omit(total)
dimnames(total)[[2]] <- c("d.y1", "d.y2", "d.y3", "l.y1", "l.y2", "l.y3", "l.x1", "l.x2", "d.y1.1", "d.y2.1", "d.y3.1", "d.x1", "d.x2" ,"d.x1.1", "d.x2.1")

y <- t(total[, 1:n])
ect <- t(total[, n + 1:(n.ect)])
x <- t(total[, -(1:(n + n.ect))])
X <- rbind(ect, x)

ols <- round(tcrossprod(y, X) %*% solve(tcrossprod(X)), 2)
dimnames(true) <- dimnames(ols)
ols
true

########################################################################################
# Bayesian estimator
iterations <- 10000
burnin <- 100

t <- dim(y)[2]
r <- 2

# Estimate the model
store <- iterations - burnin
if (store > (iterations - burnin)) {stop("Number of iterations must be larger than burn-in draws.")}

n <- dim(y)[1]
n.ect <- dim(ect)[1]
n.x <- dim(x)[1]
n.B <- n * n.x

n.alpha <- r * n
n.alpha.B <- n.alpha + n.B
n.beta <- r * n.ect

alpha.mu.prior <- matrix(0, n.alpha)
alpha.V.i.prior <- diag(1/9, n.alpha)
B.mu.prior <- matrix(0, n.B)
B.V.i.prior <- diag(1/9, n.B)
Omega.df.post <- t #+ r
Omega.V.prior <- diag(0, n)

B_filter <- matrix(0, n.B, t + 1)
Q_B <- diag(.0001, n.B)
y_tilde <- y * NA
y_tilde_0 <- y * NA
y_hat <- y * NA
y_alpha_0 <- y * NA
beta_filter <- matrix(0, n.beta, t + 1) 
beta <- matrix(stats::rnorm(n.beta), n.beta, t + 1)
beta_Q <- diag(.00001, n.beta)
Pi <- matrix(0, n * n.ect, t)
rho <- .999999
#Z.alpha <- matrix(NA, n * t, n.alpha)
#x.ect <- matrix(NA, r, t)

Omega <- diag(.0001, n)
Omega.i <- solve(Omega)

y_hat <- y

draws.B <- array(NA, c(n.B, t, iterations - burnin))
draws.Pi <- array(NA, c(n * n.ect, t, iterations - burnin))
if (r == 0) {draws.Pi[,] <- 0}
draws.Omega <- matrix(NA, n * n, iterations - burnin)
draws.LL <- matrix(NA, t, iterations - burnin)

if (iterations <= 100){
  update <- 1:iterations
} else {
  update <- round(seq(from = 0, to = iterations, length.out = 100)) # Refreshment rate of the progress bar
}
pb <- utils::txtProgressBar(width = 70, style = 3)

for (draw in 1:iterations){
  #x.ect <- rbind(matrix(NA, r, t), x)
  #for (i in 1:t) {
  #  x.ect[1:r, i] <- crossprod(matrix(beta[,i], n.ect), ect[,i])
  #}
  #Z.B <- kronecker(t(x.ect), diag(1, n))
  Z.B <- kronecker(t(x), diag(1, n))
  
  y_0 <- y * NA
  for (i in 1:t) {
    y_0[,i] <- y[,i] - Z.B[(i - 1) * n + 1:n, ] %*% B_filter[,i]
  }
  
  #B.V.i.prior <- diag(c(rep(1 / 10, n.alpha), rep(1/2, n.B)), n.alpha.B)
  B.V.i.prior <- diag(1 / 2, n.B)
  #B_0 <- posterior_normal_sur(y_0, Z.B, Omega.i, matrix(0, n.alpha.B), B.V.i.prior)
  B_0 <- posterior_normal_sur(y_0, Z.B, Omega.i, matrix(0, n.B), B.V.i.prior)
  y_filter <- y - matrix(Z.B %*% B_0, n)
  B_filter <- dk(y_filter, Z.B, Omega, Q_B, diag(1, n.B), matrix(0, n.B), Q_B)
  B <- B_filter + matrix(B_0, n.B, t +1)
  
  B_Q_res <- cbind(matrix(.0001, n.B) + rowSums((B[,-1] - B[,-(t + 1)])^2), t + 3)
  diag(Q_B) <- 1 / apply(B_Q_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
  
  
  #### Cointegration ####
  y_hat <- y * NA
  for (i in 1:t) {
    y_hat[,i] <- y[,i] - matrix(B[,i], n) %*% x[,i]
  }
  
  # Alpha
  
  
  # Beta
  
  alpha <- B[1:n.alpha,]
  Alpha <- alpha * 0
  for (i in 1:t) {
    alpha_svd_temp <- svd(matrix(alpha[,i], n))
    Alpha[,i] <- tcrossprod(alpha_svd_temp$u, alpha_svd_temp$v)
  }
  
  y_tilde <- y * NA
  Z.beta <- matrix(NA, n * t, n.beta)
  for (i in 1:t) {
    y_tilde[, i] <- y[,i] - matrix(B[-(1:n.alpha), i], n) %*% x[,i]
    Z.beta[(i - 1) * n + 1:n,] <- kronecker(matrix(Alpha[,i], n), t(ect[,i]))
    y_tilde_0[, i] <- y_tilde[, i] - Z.beta[(i - 1) * n + 1:n,] %*% beta_filter[,i]
  }
  
  beta_0 <- posterior_normal_sur(y_tilde_0, Z.beta, Omega.i, matrix(0, n.beta), diag(1 / 10 , n.beta))
  y_tilde_filter <- y_tilde - matrix(Z.beta %*% beta_0, n)
  beta_filter <- dk(y_tilde_filter, Z.beta, Omega, beta_Q, diag(1, n.beta), matrix(0, n.beta), beta_Q)
  beta <- beta_filter + matrix(beta_0, n.beta, t + 1)
  #beta <- dk(y_tilde, Z.beta, Omega, diag(1, n.beta), diag(rho, n.beta), matrix(0, n.beta), diag(1 / (1 - rho^2), n.beta))
  
  beta_Q_res <- cbind(matrix(.0001, n.beta) + rowSums((beta[,-1] - beta[,-(t + 1)])^2), t + 3)
  diag(beta_Q) <- 1 / apply(beta_Q_res, 1, function(x) {return(stats::rgamma(n = 1, shape = x[2] / 2, rate = x[1] / 2))})
  
  y_dot <- y * 0
  for (i in 1:t) {
    svd_temp <- svd(matrix(beta[,i], n.ect))
    alpha_temp <- matrix(Alpha[,i], n) %*% svd_temp$v %*% tcrossprod(diag(svd_temp$d, r), svd_temp$v)
    #beta_temp <- tcrossprod(svd_temp$u, svd_temp$v)
    beta[,i] <- tcrossprod(svd_temp$u, svd_temp$v)
    
    Pi[, i] <- tcrossprod(alpha_temp, matrix(beta[,i], n.ect))
    #Pi[, i] <- tcrossprod(matrix(alpha[,i], n), matrix(beta[,i], n.ect))
    y_dot[, i] <- y_tilde[,i] - matrix(Pi[,i], n) %*% ect[,i]
  }

  # Error covariance matrix
  Omega.i <- stats::rWishart(1, Omega.df.post, solve(tcrossprod(y_dot)))[,,1]
  Omega <- solve(Omega.i)
  
  # Likelihood calculations
  if (draw > burnin){
    pos.draw <- draw - burnin
    draws.B[,, pos.draw] <- B[-(1:n.alpha),-(t + 1)]
    draws.Pi[, ,pos.draw] <- Pi[,-(t + 1)]
    draws.Omega[, pos.draw] <- Omega
    draws.LL[, pos.draw] <- getLL(y_dot, Omega, Omega.i)
    
    if (TRUE) {
      if (draw %in% update & pos.draw > 10) {
        graphics::par(mfcol = c(2, 2))
        #if (tvp[1]) {
          stats::plot.ts(cbind(0, t(apply(draws.B[1, , (ceiling(pos.draw) * .5):pos.draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                         plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Short-run parameter") 
        #}
        #if ("const" %in% dimnames(x)[[1]]) {
          stats::plot.ts(cbind(0, Pi_true[5, pos[-(1:2)]], t(apply(draws.Pi[5, , (ceiling(pos.draw) * .5):pos.draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                         plot.type = "single", col = c("black", "blue", "red", "black", "red"), ylab = "Short-run parameter\nDeterministic") 
        #}
        #if ("const" %in% dimnames(ect)[[1]]) {
          stats::plot.ts(cbind(0, Pi_true[15, pos[-(1:2)]],t(apply(draws.Pi[15, , (ceiling(pos.draw) * .5):pos.draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
                         plot.type = "single", col = c("black", "blue", "red", "black", "red"), ylab = "Short-run parameter\nDeterministic") 
        #}
        #if (tvp[2]) {
        #  stats::plot.ts(cbind(0, t(apply(draws_Pi[8, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #                 plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Cointegration parameter") 
        #}
        #if (sv | tvp[3]) {
        #  stats::plot.ts(log(t(apply(draws_Omega[1, , (ceiling(pos_draw) * .3):pos_draw], 1, stats::quantile, probs = c(.05, .5, .95)))),
        #                 plot.type = "single", col = c("black", "red", "black", "red"), ylab = "Log-Volatility") 
        #}
        graphics::par(mfcol = c(1, 1))
      }
    }
  }
  # Update progress bar
  if (is.element(draw,update)) {utils::setTxtProgressBar(pb, draw/iterations)}
}
close(pb)


bayes <- round(cbind(matrix(rowMeans(draws.Pi), n), matrix(rowMeans(draws.B[,1,]), n)), 2)
dimnames(bayes) <- dimnames(ols)

matrix(c(bayes, ols, true), ncol = 3)
true
round(matrix(rowMeans(draws.B.shrinkage), n), 2)

test.stat <- c("AIC" = 2 * sum(log(rowMeans(draws.LL))) - (n.alpha.B) * 2,
               "BIC" = 2 * sum(log(rowMeans(draws.LL))) - (n.alpha.B) * log(t),
               "HQ" = 2 * sum(log(rowMeans(draws.LL))) - (n.alpha.B) * 2 * log(log(t)))
test.stat
