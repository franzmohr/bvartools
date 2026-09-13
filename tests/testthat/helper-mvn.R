# Log-likelihood of multivariate normal residuals, one row per period, under a
# precision matrix. Written out without the package, so that tests can check
# what the package computes against it.
mvn_loglik <- function(u, precision) {
  u <- as.matrix(u)
  k <- ncol(u)
  n <- nrow(u)
  -0.5 * n * k * log(2 * pi) +
    0.5 * n * as.numeric(determinant(precision, logarithm = TRUE)[["modulus"]]) -
    0.5 * sum((u %*% precision) * u)
}
