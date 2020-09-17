#' Summarising Bayesian VAR Models
#'
#' summary method for class \code{"bvarlist"}.
#'
#' @param object an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{draw_posterior}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The log-likelihood for the calculation of the information criteria is obtained by
#' \deqn{LL = \frac{1}{R} \sum_{i = 1}^{R} \left( \sum_{t = 1}^{T} -\frac{K}{2} \ln 2\pi - \frac{1}{2} \ln |\Sigma_t^{(i)}| -\frac{1}{2} (u_t^{{(i)}\prime} (\Sigma_t^{(i)})^{-1} u_t^{(i)} \right)},
#' where \eqn{u_t = y_t - \mu_t}. The Akaike, Bayesian and Hannan–Quinn (HQ) information criteria are calculated as
#' \deqn{AIC = 2 (Kp + Ms + N) - 2 LL},
#' \deqn{BIC = (Kp + Ms + N) ln(T) - 2 LL} and 
#' \deqn{HQ = 2  (Kp + Ms + N) ln(ln(T)) - 2 LL}, respectively.
#'
#' @return \code{summary.bvarlist} returns a table of class \code{"summary.bvarlist"}.
#'
#' @export
summary.bvarlist <- function(object, ...){
  
  n_models <- length(object)
  
  teststats <- data.frame(p = rep(NA, n_models),
                          s = rep(NA, n_models),
                          r = rep(NA, n_models),
                          LL = rep(NA, n_models),
                          AIC = rep(NA, n_models),
                          BIC = rep(NA, n_models),
                          HQ = rep(NA, n_models),
                          stringsAsFactors = FALSE)
  
  endo_vars <- NULL
  ect_vars <- NULL
  exog_vars <- NULL
  
  for (i in 1:n_models) {
    
    tt <- nrow(object[[i]][["y"]])
    k <- object[[i]]$specifications$dims["K"]
    p <- object[[i]]$specifications$lags["p"]
    endo_vars <- c(endo_vars, dimnames(object[[i]][["y"]])[[2]])
    
    teststats[i, "p"] <- p
    teststats[i, "s"] <- object[[i]]$specifications$lags["s"]
    
    temp_pars <- NULL
    x <- NULL
    
    if ("bvar" %in% class(object[[i]])) {
      
      type <- "VAR"
      x <- t(object[[i]][["x"]])
      tot_pars <- NCOL(object[[i]][["x"]])
      exog_vars <- c(exog_vars, dimnames(object[[i]][["x"]])[[2]])
      
      vars <- c("A", "B", "C")
      for (j in vars) {
        if (!is.null(object[[i]][[j]])) {
          temp_pars <- cbind(temp_pars, object[[i]][[j]])
        }
      }
    }
    
    if ("bvec" %in% class(object[[i]])) {
      
      type <- "VEC"
      if (is.null(object[[i]][["r"]])) {
        if (is.null(object[[i]][["alpha"]])) {
          r <- 0
          teststats[i, "r"] <- 0
        } else {
          r <- NCOL(object[[i]][["alpha"]]) / k
          teststats[i, "r"] <- r
        } 
      } else {
        r <- object[[i]][["r"]]
      }
      tot_pars <- r
      
      vars <- c("Pi", "Pi_x", "Pi_d", "Gamma", "Upsilon", "C")
      for (j in vars) {
        if (!is.null(object[[i]][[j]])) {
          temp_pars <- cbind(temp_pars, object[[i]][[j]])
          if (j == "Pi") {
            x <- cbind(x, object[[i]][["w"]])
            ect_vars <- c(ect_vars, dimnames(object[[i]][["w"]])[[2]])
          }
          if (j == "Pi_x") {
            x <- cbind(x, object[[i]][["w_x"]])
            ect_vars <- c(ect_vars, dimnames(object[[i]][["w_x"]])[[2]])
          }
          if (j == "Pi_d") {
            x <- cbind(x, object[[i]][["w_d"]])
            ect_vars <- c(ect_vars, dimnames(object[[i]][["w_d"]])[[2]])
          }
          if (j == "Gamma") {
            x <- cbind(x, object[[i]][["x"]])
            exog_vars <- c(exog_vars, dimnames(object[[i]][["x"]])[[2]])
          }
          if (j == "Upsilon_x") {
            x <- cbind(x, object[[i]][["x_x"]])
            exog_vars <- c(exog_vars, dimnames(object[[i]][["x_x"]])[[2]])
          }
          if (j == "C") {
            x <- cbind(x, object[[i]][["x_d"]])
            exog_vars <- c(exog_vars, dimnames(object[[i]][["x_d"]])[[2]])
          }
        }        
      }
      
      tot_pars <- tot_pars + NCOL(x)
      x <- t(x)
    }
    
    draws <- nrow(temp_pars)
    LL <- matrix(NA, tt, draws) # Get LogLik
    
    for (j in 1:draws) {
      # Residuals
      if (type == "VAR") {
        u <- t(object[[i]][["y"]]) - matrix(temp_pars[j, ], k) %*% x
      }
      if (type == "VEC") {
        u <- t(object[[i]][["y"]]) - matrix(temp_pars[j, ], k) %*% x
      }
      # Sigma
      sigma <- matrix(object[[i]]$Sigma[j,], k)
      # LogLik
      LL[, j] <- loglik_normal(u, sigma)
    }
    
    ll <- colSums(LL)
    aic <- 2 * tot_pars - 2 * ll
    bic <- log(tt) * tot_pars - 2 * ll
    hq <- 2 * log(log(tt)) * tot_pars - 2 * ll
    teststats[i, "LL"] <- mean(ll)
    teststats[i, "AIC"] <- mean(aic)
    teststats[i, "BIC"] <- mean(bic)
    teststats[i, "HQ"] <- mean(hq)
  }
  
  endo_vars <- unique(endo_vars)
  if (!is.null(ect_vars)) {
    ect_vars <- unique(ect_vars)
  }
  if (!is.null(exog_vars)) {
    exog_vars <- unique(exog_vars) 
  }

  # Omit unnecessary columns
  teststats <- teststats[, which(!apply(teststats, 2, function(x) {all(is.na(x))}))]
  
  # result <- list(model = list(endogen = endo_vars,
  #                             ect = ect_vars,
  #                             exogen = exog_vars),
  #                teststats = teststats)
  result <- teststats
  
  # class(result) <- "summary.bvarlist"
  return(result)
}