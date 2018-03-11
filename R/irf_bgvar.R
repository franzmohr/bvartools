#' Generalised Impulse Response of a GVAR Model
#' 
#' Computes the generalised impulse response (GIR) coefficients of a GVAR model.
#' 
#' @param data a list containing a solved GVAR model as produced by \code{\link{solve_gvar}}.
#' @param n.ahead an integer specifying the steps.
#' @param impulse a character vector containing the country name and variable of the impulse.
#' @param response a character vector containing the country name and variable of the response.
#' @param shock a character or numeric specifying the size of the shock. Either "sd" (default), "nsd" for negative standard deviation of the impulse variable or a numeric value.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the credible intervals.
#' @param draws a numeric specifying the number of draws from the GVAR model used to calculate the GIR.
#' @param t an integer specifying the period for which the GIR should be calculated. Only used if a GVAR with time-varying
#' parameters or stochastic volatility was estimated.
#' @param mean a logical. If \code{TRUE}, the means of the coefficient draws will be used to calculate a
#' single impulse response with no credible intervals.
#' @param cumulative a logical. If \code{TRUE}, a cumulative GIR will be calculated.
#' 
#' @return A list containing the following elements:
#' \item{IR}{a "zoo" object containing the median and credible intervals of the GIR.}
#' \item{impulse}{a character vector containing the country name and variable of the impulse.}
#' \item{response}{a character vector containing the country name and variable of the response.}
#' 
#' @export
irf_bgvar <- function(data, n.ahead = 20, impulse, response, shock = "sd", ci = .95, draws = NULL,
                      t = NULL, mean = FALSE, cumulative = FALSE){
  if (ci <= 0 | ci >= 1) {stop("Credible intervals must take a velue between 0 and 1.")}
  n <- sqrt(dim(data$GVAR$coefs$Sigma)[1])
  tvp <- dim(data$GVAR$coefs$G0.i)[2] > 1
  sv <- dim(data$GVAR$coefs$Sigma)[2] > 1
  t.max <- max(dim(data$GVAR$coefs$G0.i)[2], dim(data$GVAR$coefs$Sigma)[2])
  t.force <- FALSE
  if (is.null(t)){
    if (tvp){
      t <- t.max
    } else {
      t <- 1
    }
    if (sv){
      t.sv <- t.max
    } else {
      t.sv <- 1
    }
  } else {
    if (t > t.max | t < 1){
      stop("Specified time index is not available.")
    }
    t.force <- TRUE
    if (sv){
      t.sv <- t
    } else {
      t.sv <- 1
    }
    
    if (!tvp){
      t <- 1
    }
  }
  
  if (is.null(draws)){
    draws <- dim(data$GVAR$coefs$G0.i)[3]
  }
  
  index <- data$GVAR$specs$index
  global <- !is.na(data$GVAR$coefs$lags["global"])
  
  draws.ir <- matrix(NA, 1 + n.ahead, draws)
  if (mean){
    G <- matrix(rowMeans(data$GVAR$coefs$G[, t,]), n)
    S.u <- matrix(rowMeans(data$GVAR$coefs$Sigma[, t.sv,]), n)
    G0.i <- matrix(rowMeans(data$GVAR$coefs$G0.i[, t,]), n)
    IR <- zoo::zooreg(gir(G = G, G0.i = G0.i, Sigma = S.u, n.ahead = n.ahead, impulse = impulse, response = response,
                 index = index, shock = shock), start = 0, frequency = stats::frequency(data$GVAR$data$X))
  } else {
    pb <- utils::txtProgressBar(width = 70, style = 3)
    for (draw in 1:draws){
      G <- matrix(data$GVAR$coefs$G[, t, draw], n)
      Sigma <- matrix(data$GVAR$coefs$Sigma[, t.sv, draw], n)
      G0.i <- matrix(data$GVAR$coefs$G0.i[, t, draw], n)
      Psi <- gir(G = G, G0.i = G0.i, Sigma = Sigma, n.ahead = n.ahead, impulse = impulse, response = response,
                 index = index, shock = shock)
      
      draws.ir[, draw] <- Psi
      utils::setTxtProgressBar(pb, draw/draws)
    }
    
    ci <- 1 - ci
    ci <- apply(draws.ir, 1, stats::quantile, probs = c(.5, ci/2, 1 - ci/2))
    IR <- t(ci)
    if (cumulative){
      IR <- apply(IR, 2, cumsum) 
    }
    
    if (tvp | sv | t.force) {
      ind.start <- zoo::index(data$GVAR$data$X)[max(c(t, t.sv))]
      IR <- zoo::zooreg(IR, start = ind.start, frequency = stats::frequency(data$GVAR$data$X))
    } else {
      IR <- zoo::zooreg(IR, start = 0, frequency = stats::frequency(data$GVAR$data$X))
    }
  }
  
  result <- list(IR = IR, impulse = impulse, response = response)
  return(result)
}