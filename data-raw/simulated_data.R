t <- 200
p <- 1

tt <- t + p * 2 + 1
a11 <- rep(.5, tt) - 1:tt * .002
a22 <- rep(-.4, tt) + 1:tt * .001

y <- matrix(0, 2, tt)

set.seed(7895)
for (i in (p + 1):tt) {
  y[, i] <- matrix(c(a11[i], 0, 0, a22[i]), 2) %*% y[, i - 1] + rnorm(2, 0, 1)
}

y <- t(y[, -(1:(p + 1))])
#plot.ts(t(y))

