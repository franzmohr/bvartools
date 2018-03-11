#' Add Minnesota Priors to the Country Models
#' 
#' 
#' 
#' @param data a list as produced by \link{\code{country_models}}.
#' @param method
#' 
#' @references 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection. \emph{Journal of Applied Econometrics}, 28(2), 204--230.
#' @export
add_Minnesota_prior <- function(data, kappa = 1, lambda = 1){
  if (length(lambda) == 1) {
    lambda = c(1, 1)
  }
  
  priors <- lapply(data$Country.data, function(x) {
    if (x$specs$type == "VAR") {
      temp <- gen.varx(x)
      B <- tcrossprod(temp$y, temp$x)%*%solve(tcrossprod(temp$x))
      res <- temp$y - B%*%temp$x
      s2 <- matrix(diag(tcrossprod(res) / (ncol(temp$y) - 1)))
      s2.foreign <- matrix(diag(tcrossprod(temp$x - rowMeans(temp$x)) / (ncol(temp$x) - 1)))
      
      n <- temp$domestic["dim"]
      n.x <- temp$foreign["dim"]
      n.g <- temp$global["dim"]
      p <- temp$domestic["lag"]
      q <- temp$foreign["lag"]
      s <- temp$global["lag"]
      
      prior <- B * 0
      # Own lags
      for (i in 1:p) {
        prior[, (i-1)*n + 1:n] <- lambda[1] * tcrossprod(s2, 1/s2) / i^2
        diag(prior[, (i-1)*n + 1:n]) <- kappa / i^2
      }
      for (i in 1:(q + 1)) {
        prior[, n * p + (i-1) * n.x + 1:n.x] <- lambda[2] * tcrossprod(s2, 1 / s2.foreign[n * p + (i-1) * n.x + 1:n.x,]) / i^2
      }
      if (n.g > 0){
        for (i in 1:(s + 1)) {
          prior[, n * p + n.x * q + (i-1) * n.g + 1:n.g] <- lambda[2] * tcrossprod(s2, 1 / s2.foreign[n * p + n.x * q + (i-1) * n.g + 1:n.g,]) / i^2
        }
      }
      
      if (x$specs$deterministic.terms %in% c("II", "III", "IV", "V")) {
        prior[, temp$domestic["total"] + temp$foreign["total"] + temp$global["total"] + 1] <- 100 * s2
      }
      if (x$specs$deterministic.terms %in% c("IV", "V")) {
        prior[, temp$domestic["total"] + temp$foreign["total"] + temp$global["total"] + 2] <- 100 * s2
      }
      
      prior <- diag(c(prior))
      return(prior)
    }
  })

  for (i in 1:length(data$Country.data)) {
    data$Country.data[[i]]$priors$A$constant[[2]] <- priors[[i]]
  }
  
  return(data)
}