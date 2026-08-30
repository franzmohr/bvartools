#' Summarising Bayesian VEC Coefficients
#'
#' summary method for class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param period integer. Index of the period, for which the summary statistics should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param x an object of class 'summary.bvecmodel', usually, a result of a call to
#' \code{\link{summary.bvecmodel}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{summary.bvecmodel} returns a list of class 'summary.bvecmodel',
#' which contains the following components:
#' \item{a}{A list of various summary statistics of the posterior
#' draws of the VEC coefficients.}
#' \item{sigma}{A list of various summary statistics of the posterior
#' draws of the variance-covariance matrix.}
#' \item{model}{a list containing information on the model specification.}
#'
#' @export
summary.bvecmodel <- function(object, ci = .95, period = NULL, ...){
  
  # Number of endogenous variables
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  r <- object[["model"]][["rank"]]
  tt <- nrow(object[["data"]][["train"]][["y"]])
  use_incl <- object[["model"]][["varsel"]] %in% c("ssvs", "bvs")
  structural <- object[["model"]][["structural"]]
  sv <- object[["model"]][["error"]] %in% c("sv", "sv+covar")
  tvp <- object[["model"]][["tvp"]]
  tvp_and_covar <- tvp & object[["model"]][["error"]] %in% c("gamma", "gamma+covar")
  if (sv | tvp) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
  }
  

  y_names <- object[["model"]][["endogen"]]
  ect_names <- dimnames(object[["data"]][["train"]][["w"]])[[2]]
  x_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]]
  if (r > 0) {
    x_names <- c(ect_names, x_names)
  }
  dim_names <- list(y_names, x_names)
  
  # Non-error coefficients
  means <- NULL
  median <- NULL
  sds <- NULL
  naive_sd <- NULL
  ts_sd <- NULL
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  q_low <- NULL
  q_high <- NULL
  
  result <- NULL
  result[["model"]] <- object[["model"]]
  
  if (r > 0) {
    n_alpha <- r * k
    k_ect <- ncol(object[["data"]][["train"]][["w"]])
    draws <- nrow(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    Pi <- matrix(NA_real_, draws, k * k_ect)
    for (i in 1:draws) {
      Pi[i, ] <- tcrossprod(matrix(object[["posterior"]][["a"]][["coeffs"]][i, 1:n_alpha], k),
                            matrix(object[["posterior"]][["beta"]][["coeffs"]][i, ], ncol = r))
    }
    temp_par <- coda::mcpar(object[["posterior"]][["a"]][["coeffs"]])
    Pi <- coda::mcmc(Pi, start = temp_par[1], end = temp_par[2], thin = temp_par[3])
    
    temp <- summary(Pi, quantiles = c(ci_low, .5, ci_high))
    if ("numeric" %in% class(temp[["statistics"]])) {
      means <- cbind(means, matrix(temp[["statistics"]]["Mean"], k))
      sds <- cbind(sds, matrix(temp[["statistics"]]["SD"], k))
      naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]]["Naive SE"], k))
      ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]]["Time-series SE"], k))
      q_low <- cbind(q_low, matrix(temp[["quantiles"]][1], k))
      median <- cbind(median, matrix(temp[["quantiles"]][2], k))
      q_high <- cbind(q_high, matrix(temp[["quantiles"]][3], k)) 
    } else {
      means <- cbind(means, matrix(temp[["statistics"]][, "Mean"], k))
      sds <- cbind(sds, matrix(temp[["statistics"]][, "SD"], k))
      naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]][, "Naive SE"], k))
      ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]][, "Time-series SE"], k))
      q_low <- cbind(q_low, matrix(temp[["quantiles"]][, 1], k))
      median <- cbind(median, matrix(temp[["quantiles"]][, 2], k))
      q_high <- cbind(q_high, matrix(temp[["quantiles"]][, 3], k)) 
    }
  }
  
  if (!is.null(object[["posterior"]][["a"]][["coeffs"]])) {
    
    n_a <- ncol(object[["data"]][["train"]][["z"]])
    pos <- 1:n_a
    if (r > 0) {
      pos <- pos[-(1:n_alpha)]
    }
    if (tvp) {
      pos <- (period - 1) * n_a + pos
    }
    
    temp <- summary(object[["posterior"]][["a"]][["coeffs"]][, pos], quantiles = c(ci_low, .5, ci_high))
    
    if ("numeric" %in% class(temp[["statistics"]])) {
      means <- cbind(means, matrix(temp[["statistics"]]["Mean"], k))
      sds <- cbind(sds, matrix(temp[["statistics"]]["SD"], k))
      naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]]["Naive SE"], k))
      ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]]["Time-series SE"], k))
      q_low <- cbind(q_low, matrix(temp[["quantiles"]][1], k))
      median <- cbind(median, matrix(temp[["quantiles"]][2], k))
      q_high <- cbind(q_high, matrix(temp[["quantiles"]][3], k)) 
    } else {
      
      if (structural) {
        
        stop("update structural")
        n_non_structural <- k * (k * p + m * (s + 1) + n)
        struct_matrix <- matrix(1:(k * k), k)
        pos_values <- which(lower.tri(struct_matrix))
        pos_zero <- which(upper.tri(struct_matrix))
        pos_one <- struct_matrix[-c(pos_values, pos_zero)] + n_non_structural
        pos_values <- pos_values + n_non_structural
        pos_zero <- pos_zero + n_non_structural
        
        pos_a <- matrix(NA, k , k)
        n_structural <- k * (k - 1) / 2
        pos_a[upper.tri(pos_a)] <- 1:n_structural
        pos_a <- t(pos_a)
        pos_a <- pos_a[lower.tri(pos_a)] + n_non_structural
        
        res_statistics <- matrix(NA_real_, n_non_structural + k * k, ncol(temp[["statistics"]]))
        dimnames(res_statistics) <- list(NULL, dimnames(temp[["statistics"]])[[2]])
        res_statistics[1:n_non_structural,] <- temp[["statistics"]][1:n_non_structural,]
        res_statistics[pos_one, 1] <- 1
        res_statistics[pos_one, -1] <- 0
        res_statistics[pos_zero, ] <- 0
        res_statistics[pos_values, ] <- temp[["statistics"]][pos_a,]
        temp[["statistics"]] <- res_statistics
        
        res_quantiles <- matrix(NA_real_, n_non_structural + k * k, ncol(temp[["quantiles"]]))
        dimnames(res_quantiles) <- list(NULL, dimnames(temp[["quantiles"]])[[2]])
        res_quantiles[1:n_non_structural,] <- temp[["quantiles"]][1:n_non_structural,]
        res_quantiles[pos_one, ] <- 1
        res_quantiles[pos_zero, ] <- 0
        res_quantiles[pos_values, ] <- temp[["quantiles"]][pos_a,]
        temp[["quantiles"]] <- res_quantiles
        
      }
      
      means <- cbind(means, matrix(temp[["statistics"]][, "Mean"], k))
      sds <- cbind(sds, matrix(temp[["statistics"]][, "SD"], k))
      naive_sd <- cbind(naive_sd, matrix(temp[["statistics"]][, "Naive SE"], k))
      ts_sd <- cbind(ts_sd, matrix(temp[["statistics"]][, "Time-series SE"], k))
      q_low <- cbind(q_low, matrix(temp[["quantiles"]][, 1], k))
      median <- cbind(median, matrix(temp[["quantiles"]][, 2], k))
      q_high <- cbind(q_high, matrix(temp[["quantiles"]][, 3], k)) 
    }
    
    if (use_incl) {
      incl <- colMeans(object[["posterior"]][["a"]][["lambda"]])
      if (structural) {
        res <- rep(NA, n_non_structural + k * k)
        res[1:n_non_structural] <- incl[1:n_non_structural]
        res[pos_one] <- 1
        res[pos_zero] <- 0
        res[pos_values] <- incl[pos_a]
        incl <- res
      }
      incl <- matrix(incl, k) 
    }
    
    if (!is.null(means)) {
      dimnames(means) <- dim_names
      dimnames(sds) <- dim_names
      dimnames(naive_sd) <- dim_names
      dimnames(ts_sd) <- dim_names
      dimnames(q_low) <- dim_names
      dimnames(median) <- dim_names
      dimnames(q_high) <- dim_names
      if (use_incl) {
        dimnames(incl) <- dim_names
      }
    }
    
    result[["a"]] <- list(means = means,
                          median = median,
                          sd = sds,
                          naivesd = naive_sd,
                          tssd = ts_sd,
                          q_lower = q_low,
                          q_upper = q_high)
    
    if (use_incl) {
      result[["a"]][["lambda"]] = incl
    }  
  }
  
  
  dim_names <- list(dim_names[[1]], dim_names[[1]])
  
  # Error coefficients
  if (!is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    
    if (sv | tvp_and_covar) {
      temp <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]][, (tt - 1) * k * k + 1:(k * k)]
    } else {
      temp <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]]
    }
    
    if (k == 1) {
      temp <- matrix(1 / temp)
    } else {
      temp <- t(apply(temp, 1, function(x, k) {solve(matrix(x, k))}, k = k))
    }
    
    temp_tsp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    temp <- coda::mcmc(temp, start = temp_tsp[1], end = temp_tsp[2], thin = temp_tsp[3])
    rm(temp_tsp)
    # Summary of an mcmc object
    temp <- summary(temp, quantiles = c(ci_low, .5, ci_high))
    
    if (k == 1) {
      means <- matrix(temp[["statistics"]]["Mean"], k)
      sds <- matrix(temp[["statistics"]]["SD"], k)
      naive_sd <- matrix(temp[["statistics"]]["Naive SE"], k)
      ts_sd <- matrix(temp[["statistics"]]["Time-series SE"], k)
      q_low <- matrix(temp[["quantiles"]][1], k)
      median <- matrix(temp[["quantiles"]][2], k)
      q_high <- matrix(temp[["quantiles"]][3], k)
    } else {
      means <- matrix(temp[["statistics"]][, "Mean"], k)
      sds <- matrix(temp[["statistics"]][, "SD"], k)
      naive_sd <- matrix(temp[["statistics"]][, "Naive SE"], k)
      ts_sd <- matrix(temp[["statistics"]][, "Time-series SE"], k)
      q_low <- matrix(temp[["quantiles"]][, 1], k)
      median <- matrix(temp[["quantiles"]][, 2], k)
      q_high <- matrix(temp[["quantiles"]][, 3], k)
    }
    
    if (!is.null(object[["posterior"]][["psi"]][["lambda"]])) {
      incl <- matrix(colMeans(object[["posterior"]][["psi"]][["lambda"]]), k)
    }
    
    
    dimnames(means) <- dim_names
    dimnames(sds) <- dim_names
    dimnames(naive_sd) <- dim_names
    dimnames(ts_sd) <- dim_names
    dimnames(q_low) <- dim_names
    dimnames(median) <- dim_names
    dimnames(q_high) <- dim_names
    
    result[["sigma"]] <- list(means = means,
                              median = median,
                              sd = sds,
                              naivesd = naive_sd,
                              tssd = ts_sd,
                              q_lower = q_low,
                              q_upper = q_high)
    
    if (!is.null(object[["posterior"]][["psi"]][["lambda"]])) {
      dimnames(incl) <- dim_names
      result[["sigma"]][["lambda"]] = incl
    }
  }
  
  result[["model"]][["ci"]] <- paste(c(ci_low, ci_high) * 100, "%", sep = "")
  result[["model"]][["period"]] <- period
  
  class(result) <- list("summary.bvecmodel", "list")
  return(result)
}
